# XL / XXL pending windows — pre-registered by the weekly, queued by the daily

**Why this file exists (operator directive 2026-09-13).** XL and XXL windows
are cron scripts with no Claude session — the only compute in the project
that costs no tokens — but only the weekly review may commission one, the
queue files hold one item each, and the weekly runs at most twice a
week (once a week during the 2026-09-13 … 09-26 wind-down). Left as it was, at most one of the three budgeted XL
windows and, some weeks, no XXL window would run. This file splits the role:
the **weekly decides** (it pre-registers entries here — tier, chunk and step,
the exact command, the measured price, the readout and the decision rule)
and the **daily review is the clerk** (daily-review.md step 6b: it copies
*every* READY entry of a tier into its own file in
`docs/testing/<tier>-queue.d/` verbatim, up to the tier's budget — a FIFO
the launcher drains one window per night, 2026-09-15). An interactive
operator session may do either.

**Backlog floor and daily licence (operator directive 2026-09-15).** The
daily review keeps **≥ 4 `xl` and ≥ 1 `xxl` entries ahead** (READY or
QUEUED for a future night), reports the count on the dashboard, and when
short may write entries itself in two classes only — **priced-family
variants** of an already-measured case (same script and mesh; frequency,
degree, drive/port set, rank count or one env knob varied) and **cost
probes** (first window on an unpriced case, readout = the price). Each such
entry names its licence class. Everything else is the weekly's.

**Budget, as the clerk reads it.** `xl`: fewer than six rows with a
non-zero `Elapsed` in the trailing 7 days of `docs/testing/xl-ledger.md`
(nightly Sun–Fri, operator directive 2026-09-13)
(killed and 0-second rows are uncharged, §5.1). `xxl`: no such row in the
trailing 7 days of `docs/testing/xxl-ledger.md`. The next run date is the
next 02:00 for the tier (xl Sun–Fri, xxl Saturday) after the budget opens.

**Status words.** `READY` — prerequisite landed, command final. `PENDING
PREREQUISITE` — names the item that must land first (an ordinary §9 item).
`QUEUED <date> for <run date>` — copied into `<tier>-queue.d/`. `RUN <log>` —
the launcher consumed the queue file; the review that finds the ledger row marks
it and the weekly interprets it. Entries are never edited after `QUEUED`
except to add the run status; a changed command is a new entry.

## Entries

### 1. `xxl` — `WF-7` step 0b: the F-human rung at degree 2, one drive, 64 MHz

**Status:** RUN `20260919T070008Z_WF-7-step0b.log` — both legs rc 0,
Status 0, 673 s, `memory.peak` 106.1 GiB; degree 2 is 3 256 418 unknowns,
347.56 s, 107.9 GiB summed; `|S₂ − S₁| / |S₁|` = 5.50 % from the two printed
`S_driven` values (ledger row 2026-09-19; marked by the 2026-09-19 03:00
review; the weekly of 2026-09-19 reads it against the decision rule below —
branch (c) is excluded, and the F-small 64 MHz moves it is compared with are
6.51 / 2.60 / 4.11 %, entry 4's "Why"). *(Was QUEUED 2026-09-13 (operator,
interactive) for Saturday 2026-09-19 02:00 — moved to
`docs/testing/xxl-queue.d/20260919-WF-7-step0b.env` 2026-09-15 when the queue
became a directory; consumed by the launcher.)*

**Prerequisite:** `scripts/probes/wf7_step0_f_human_cost.py` gained the
`FEM_EM_WF7_DEGREE` knob (this commit; unset or `1` is byte-identical to
step 0). Landed.

**What runs:** the step-0 probe twice in one window against
`fem-em-solver-xxl` at `-n 16`: degree 1 first (the control — step 0's
solve at `-n 8` read 607 039 unknowns, 37 s, 10.93 GiB summed; the same
mesh at `-n 16` re-reads `S_driven`), then degree 2. Durable capture, both
`rc`s recorded, `memory.peak` printed, `[capture] rc=` last.

**Price, predicted (not measured at this order and scale — that is the
point):** degree 2 on the 507 266-cell rung is ≈ 3.2 M unknowns. From
`ANS-4` step 2d's measured point (3.79 M unknowns, 290.2 GiB peak,
≈ 1 400 s per drive before factor reuse, `-n 16`) scaled at memory ∝ N^1.35
and time ∝ N^1.8: **≈ 235 GiB, ≈ 17 min**; the probe prints its own
brackets 150–350 GiB and 10–60 min and says INSIDE / OUTSIDE. Timeouts
1 500 s + 27 000 s inside the 28 800 s window; 754 GiB limit.

**Readout (record, 🧪 — nothing asserted beyond the imported cell band):**
`S_driven` at degree 1 and degree 2 on the same F-human mesh, and
`|S₂ − S₁| / |S₁|` — **the human-scale order sensitivity, which nobody has
measured**; plus the degree-2 price (unknowns, solve time, `ru_maxrss` per
rank, summed, `memory.peak`), banked in `xxl-ledger.md`.

**Decision rule, pre-registered for the 2026-09-19 weekly:** compare the
F-human degree-1 → 2 move with F-small's own move at 64 MHz (`ANS-4` step
2b, our figures). (a) F-human's move is in the same class or larger ⇒ a
human-scale S quoted at degree 1 carries that error, degree 2 is the
production order at human scale, and this window is its price. (b) Much
smaller (an order of magnitude) ⇒ degree 1 suffices at human scale for S,
and the human-scale sweep stays an ordinary heavy-tier job. (c) The
degree-2 run does not finish or exceeds the box ⇒ the price is the finding;
`TH-16` (symmetry planes) moves up. Either way no band moves and no
absolute claim is made.

### 2. `xl` — `ANS-4` step 3: the 64 MHz order-matched rung (weekly 2026-09-13, §10)

**Status:** RUN `20260916T070008Z_ANS-4-step3.log` — 16 passed, Status 0,
1658 s, `memory.peak` 280.8 GiB (ledger row 2026-09-16; marked by the
2026-09-16 03:00 review; the weekly of 2026-09-19 interprets it against the
decision rule below). *(Was QUEUED 2026-09-15 (operator, interactive) for
Wednesday 2026-09-16 02:00 — `docs/testing/xl-queue.d/20260916-ANS-4-step3.env`,
consumed by the launcher. Before that PENDING PREREQUISITE — `ANS-4` step 3a. 3a landed on `main`
2026-09-14, 04:30 slot. The knob reaches the solve, and the flag-off control
reproduced 2a″ at its `-n 8` record width, not `-n 2` (`OPS-41`). Marking
this entry READY is the next daily review's job.)*

**Prerequisite (an ordinary §9 item, standard tier, tests only):**
`tests/validation/test_ans4_resolution_ladder.py` is hard-wired to
`FREQUENCY_128_HZ` (two sites). Step 3a adds an env knob
`FEM_EM_ANS4_FREQUENCY_HZ` (default 128e6; unset is bit-identical) and
proves it with the flag-off control — the `0.015:1` rung reproduces step
2a″'s recorded digits — through the harness at `-n 2`. No band, no record.

**Budget:** open — three charged rows in the trailing 7 days (09-09 ×2,
09-10) against six. The gate is 3a: the first daily review after 3a lands
on `main` queues this for the next 02:00 (3a queued Monday 09-14 and
landed that morning ⇒ the Wednesday 09-16 03:00 review queues it for
**Thursday 2026-09-17 02:00**; Tuesday has no review).

**Command (final once 3a lands; the knob name above is the contract):**

```
XL_CHUNK="ANS-4-step3"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=64e6 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3-raw.log; cat /workspace/logs/ans4-step3-raw.log; exit \$rc'"
```

**Price, measured:** step 2d ran this mesh and order at 128 MHz — 592 744
cells, 3.79 M unknowns, 290.2 GiB peak, 5 653 s for four drives before
factor reuse. With one factorisation and three back-substitutions
(`PORT-19`) the window is *predicted* 2 000–3 000 s; memory the same.

**Readout and decision rule:** as pre-registered by the 2026-09-13 weekly
(§10): public — the three C4 classes at 64 MHz at degree 2 against the
degree-1 record, every imported `PORT-11` gate on the rung; private
(`docs/private/`) — the miss against the AED First Order column. A miss in
the same class as the 128 MHz order-matched residual ⇒ the 64 MHz AGREE
stands on evidence and its "by mechanism" qualifier is dropped; a
materially larger miss ⇒ a frequency-dependent feed-model systematic
re-opens as a known-issues entry and a `PORT` systematics chunk.

### 3. `xl` — `ANS-4` step 3b: the 128 MHz degree-2 rung on the C4-congruent cut (daily licence: priced-family variant, 2026-09-16 review)

**Status:** RUN `20260917T070008Z_ANS-4-step3b.log` — 15 passed / 1
skipped (the `-n 16` record-width skip this entry predicted), Status 0,
1539 s, `memory.peak` 273.1 GiB, 592 550 cells at h = 0.005 under the cut
(ledger row 2026-09-17; marked by the 2026-09-18 03:00 review; the weekly
of 2026-09-19 interprets it against the decision rule below). *(Was QUEUED
2026-09-16 for Thursday 2026-09-17 02:00 —
`docs/testing/xl-queue.d/20260917-ANS-4-step3b.env`, consumed by the
launcher.)*

**Licence class:** priced-family variant of the `ANS-4-step3` (2026-09-16)
and `ANS-4-step2d` (2026-09-10) ledger rows — the same module, the same
`RUNGSPEC="0.015:1 0.005:2"`, the same generator; **one named env knob
varied**: `FEM_EM_ANS4_STEP2_C4_CONGRUENT=1` (`GEO-32`'s congruent
port-sheet cut, off in both priced rows). The frequency is the module
default (128 MHz, knob unset). The cut changes the mesh by construction
(116 118 vs 116 085 cells at the ×1 rung, +0.03 %; the h = 0.005 count
under the cut is unmeasured and is part of the readout).

**Why (the §10 question it answers):** the 128 MHz AGREE was banked on
step 2d's degree-2 rung with the cut *off*, while the degree-1 record it is
compared against (2a″) has the cut *on*, and `GEO-32` showed the cut
carries the whole degree-1 read-back spread. This window gives the
order-matched 4×4 on the same cut as the degree-1 record — apples to
apples — and shows whether the degree-2 class spreads (0.1985 / 0.1715 /
0.1946 % at 64 MHz, 3× the degree-1 rung's) are the cut too.

**Command:**

```
XL_CHUNK="ANS-4-step3b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3b-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3b-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3b-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3b-raw.log; cat /workspace/logs/ans4-step3b-raw.log; exit \$rc'"
```

**Price, inherited from the family (scaling stated):** `ANS-4-step3` ran
this `RUNGSPEC` at `-n 16` in **1658 s** with `memory.peak` **280.8 GiB**
(592 744 cells, 3.79 M unknowns, four drives 1476.9 s); step 2d's same
rung read 290.2 GiB. A direct solver's cost is set by the mesh and order,
not the frequency, and the cut moves the cell count by ≈0.03 %, so the
prediction is **1 500–2 500 s, 270–300 GiB**; if the congruent cut at
h = 0.005 meshes materially differently, the cell count printed at
`[ANS-4 step2d] rung h=0.005 degree 2 (c4_congruent_sheets=on)` is the
first reading. Timeout 14 400 s inside the 4 h window; 512 GiB limit.

**Readout (record, nothing asserted beyond the module's imported gates):**
public — the degree-2 4×4 at 128 MHz on the congruent cut, its class
spreads beside step 2d's (cut off), and the degree 1 → 2 move on the same
cut (the ×1 rung in this window against the 2a″ record, which the module
prints); private (`docs/private/`) — the miss against the AED First Order
column beside step 2d's. **Decision rule for the 09-19 weekly:** classes
within their own step-2d spread of step 2d's ⇒ the cut is not a degree-2
systematic and the 128 MHz AGREE stands as written; a class move larger
than the spread ⇒ the cut is a named systematic of the order-matched rung
and `GEO-32`'s 🧪 row gains the degree-2 reading. No band moves. The
module's 128 MHz record assert (`test_the_frequency_knob_reaches_the_solve`)
is skipped at `-n 16` by design (`OPS-41`) — the ×1 digits are printed.

### 4. `xl` — `ANS-4` step 3c: the 10 MHz degree-2 rung (daily licence: priced-family variant, 2026-09-16 review)

**Status:** RUN `20260918T070008Z_ANS-4-step3c.log` — 16 passed, Status 0,
1603 s, `memory.peak` 284.8 GiB (ledger row 2026-09-18; marked by the
2026-09-18 03:00 review; the weekly interprets it against the decision
rule below — and reads the self-class spread, 0.4451 % against the 0.5 %
band, with entry 8). *(Was QUEUED 2026-09-16 for Friday 2026-09-18 02:00 —
`docs/testing/xl-queue.d/20260918-ANS-4-step3c.env`, consumed by the
launcher.)*

**Licence class:** priced-family variant of `ANS-4-step3` — same module,
same `RUNGSPEC`, same mesh (cut off, the `GEO-19` record mesh), **frequency
varied only**: `FEM_EM_ANS4_FREQUENCY_HZ=10e6`.

**Why:** the `ANS-4` diagnosis began as a *frequency trend* — AED agreed
with our degree-1 classes at 10 MHz and disagreed at the Larmor rungs. The
degree 1 → 2 move is now measured at 128 MHz (6.09 / 5.38 / 6.70 %) and
64 MHz (6.51 / 2.60 / 4.11 %). The 10 MHz point says whether the order
sensitivity itself is frequency-dependent the way the miss was: a small
move at 10 MHz means the 10 MHz AGREE was non-accidental at degree 1 and
the Larmor misses were order; a move of the same size means order matters
everywhere and the 10 MHz agreement was at a degree-1 point that happened
to sit inside AED's band.

**Command:**

```
XL_CHUNK="ANS-4-step3c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=10e6 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3c-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3c-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3c-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3c-raw.log; cat /workspace/logs/ans4-step3c-raw.log; exit \$rc'"
```

**Price, inherited:** `ANS-4-step3` — **1658 s, 280.8 GiB** at `-n 16` on
this exact mesh; frequency does not change the factorisation's size, so
the prediction is **1 500–2 000 s, 270–290 GiB**.

**Readout (record):** public — the degree-2 4×4 at 10 MHz, the imported
`PORT-9` gates on the rung, and the degree 1 → 2 move on the driven column
(the ×1 rung's 10 MHz 4×4 is `PORT-9` leg (d)'s, reproduced by `ANS-4`
run1 to 1.2e-10); private — the miss against the AED 10 MHz columns beside
the degree-1 one. **Decision rule for the weekly:** a 10 MHz move an order
of magnitude below the Larmor moves ⇒ the frequency trend is an
order-sensitivity trend, and the 10 MHz AGREE stands at either order; a
move in the same class ⇒ the degree-1 10 MHz AGREE is re-read as
coincidental and `ANS-4`'s row says so. The module's off-128 MHz knob
assert (every class moves > 1e-2 from the 128 MHz record) is the only
frequency-dependent gate and holds trivially at 10 MHz.

### 5. `xl` — `ANS-4` step 3d: the 64 MHz degree-2 rung on the C4-congruent cut (daily licence: priced-family variant, 2026-09-16 review)

**Status:** RUN `20260920T070008Z_ANS-4-step3d.log` — 16 passed, 0 skipped,
Status 0, 1619 s, `memory.peak` 281.1 GiB; class spreads on the cut
0.0132 / 0.0126 / 0.0081 %, two-rung move 6.50 / 2.61 / 4.12 % beside step
3's 6.51 / 2.60 / 4.11 % (ledger row, 2026-09-20 review). By the public rule
below the cut is not a degree-2 systematic at 64 MHz; the 09-26 weekly
rules and reads the private column. *(Was QUEUED 2026-09-16 for Sunday
2026-09-20 02:00 — `docs/testing/xl-queue.d/20260920-ANS-4-step3d.env`,
consumed by the launcher.)*

**Licence class:** priced-family variant of `ANS-4-step3` — same module,
same `RUNGSPEC`, same frequency (64 MHz); **one named env knob varied**:
`FEM_EM_ANS4_STEP2_C4_CONGRUENT=1`. The 64 MHz counterpart of entry 3, so
the weekly reads the cut's effect at both Larmor frequencies.

**Command:**

```
XL_CHUNK="ANS-4-step3d"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=64e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3d-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3d-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3d-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3d-raw.log; cat /workspace/logs/ans4-step3d-raw.log; exit \$rc'"
```

**Price, inherited:** as entry 3 — **1 500–2 500 s, 270–300 GiB** at
`-n 16` (the congruent-cut mesh at h = 0.005 will have been measured by
entry 3 the night before; if entry 3 found it materially larger, the
launcher's budget check is unaffected but the weekly should re-read this
prediction).

**Readout and decision rule:** as entry 3, at 64 MHz, beside
`ANS-4-step3`'s 4×4 (cut off): classes within step 3's own spread ⇒ the
cut is not a degree-2 systematic at 64 MHz either; a larger move ⇒ named
systematic, `GEO-32` row annotated. Private: the miss against AED's 64 MHz
First Order column beside step 3's. No band moves.

### 6. `xl` — `WF-7` step 0c: the F-human degree-1 solve at the full 32-port set (daily licence: cost probe, 2026-09-16 review)

**Status:** RUN `20260921T070008Z_WF-7-step0c.log` — `[capture] rc=0`,
Status 0, 228 s, `memory.peak` 11.96 GiB; one factorisation (24.73 s) + 31
held back-substitutions (0.40–0.57 s), cumulative 3.74 min — **below** the
5–45 min bracket, memory INSIDE; full 32×32 reciprocity 1.084e-14 asserted
green. Ledger row filled by the 2026-09-21 03:00 review. **Two corrections to
this entry's own text (log-pathologist, 2026-09-21):** the Readout's "nothing
asserted beyond the imported cell band" is wrong — the probe also asserts the
32×32 reciprocity and the drive-2 ≠ drive-1 floor (`:394–428`); and the
refinement paragraph called its own ≈ 3.5 min extrapolation "well inside"
a bracket whose floor is 5 min. Neither moves the decision rule: the set fits
one window with two orders of magnitude to spare, so chain step H3 is written
at the **heavy** tier (§9, 2026-09-21). *(Was QUEUED 2026-09-18 for Monday
2026-09-21 02:00 —
`docs/testing/xl-queue.d/20260921-WF-7-step0c.env`, the first `xl` night
after entry 5's 09-20; Saturday 09-19 is the `xxl` night; READY before that —
2026-09-16 09:38Z, §9 item 1 landed on `main` as `a267c5a`.)* The
probe's `FEM_EM_WF7_PORTS` knob is in (`scripts/probes/wf7_step0_f_human_cost.py`),
proved at heavy tier by the flag-off control (knob unset reproduces step 0's
printed `S_driven` `0.407423+0.344417j`, `20260916T093301Z_WF-7-step0c.log:10428`)
and the two-drive reciprocity assert (2×2 ratio 9.767e-16 ≤ the imported 1e-3,
drive 2 differs by 1.407e-03, `20260916T093548Z_WF-7-step0c.log:10433,10438`).
The command below is unchanged and runs as written. The next review queues it
for the following `xl` night.

**Price refinement from the two-drive window (2026-09-16, `-n 8`):** drive 1
factorises in 27.71 s, drive 2 back-substitutes the held factor in **0.59 s**
(`solve_kind` `held`), and summed `ru_maxrss` moved 10.892 → 11.082 GiB over
the second column. Extrapolated to 32 drives at `-n 8`: ≈ 2 min mesh + 28 s +
31 × 0.6 s ≈ 3.5 min and ≈ 17 GiB — well inside the 5–45 min / 10–40 GiB
brackets below, which stand unchanged (the XL window runs at `-n 16`, where
neither has been measured).

**Licence class:** cost probe — the first XL window on an unpriced case
(§5.1: "the first XL window on it is a cost probe whose readout says so").
`WF-7` step 0 priced *one* degree-1 drive on this mesh at `-n 8` on the
ordinary service (37 s solve, 123 s mesh, 10.93 GiB summed,
`20260913T190102Z_WF-7-step0.log:10419–10427`); the full 32-column set at
`-n 16` on the XL service has never run. The weekly posed the question
(daily-review.md step 6b.4: "the F-human degree-1 solve at the full port
set").

**What runs:** the step-0 probe once with every ring port driven in turn
(`FEM_EM_WF7_PORTS=all`, degree 1, 64 MHz, `-n 16`), each column's fields
dropped after its `S` column is read; the 32×32 assembled and its
reciprocity ratio, `σ_max` and class spreads printed beside `PORT-13`'s
imported bands.

**Command (final once the knob lands; the knob name is the contract):**

```
XL_CHUNK="WF-7-step0c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/wf7-step0c-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF7_PORTS=all timeout -k 60 7200 mpiexec -n 16 python3 scripts/probes/wf7_step0_f_human_cost.py >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, predicted (unmeasured — that is the point):** from step 0's
37 s per drive at `-n 8` without factor reuse, 32 drives ≈ 20 min plus
the 2 min mesh; at `-n 16` and with `reuse_factorization` after the first
drive the solve part could fall to a few minutes. **Brackets: 5–45 min,
10–40 GiB summed `ru_maxrss`**; the probe prints INSIDE / OUTSIDE per
drive and cumulatively, so a killed window still prices how far it got.
Timeout 7 200 s inside the 4 h window.

**Readout (the price, and records only):** unknowns, per-drive solve time
with and without reuse, cumulative wall clock, `ru_maxrss` per rank /
summed, `memory.peak`; the 32×32's `_reciprocity_ratio`, `σ_max` and class
spreads **printed** beside `RECIPROCITY_BAND` / `COLUMN_PASSIVITY_CEILING`
/ `OPPOSITE_SPREAD_BAND` — nothing asserted beyond the imported cell band.
**Decision rule for the weekly:** the price dates the human-scale S sweep
in §10 (a full set inside one `xl` window ⇒ the F-human degree-1 32×32 is
an ordinary XL job; outside ⇒ `TH-16` symmetry planes move up); a
reciprocity or passivity miss on the human-scale mesh is a finding for a
`PORT` chunk, not a band change.

### 7. `xl` — `ANS-4` step 3e: the 64 MHz degree-2 *h*-ladder (daily licence: priced-family variant, 2026-09-18 review)

**Status:** RUN `20260922T070007Z_ANS-4-step3e.log` — 15 passed / 2 failed,
Status 1, 2206 s, 277.44 GiB; both reds non-physics (the `OPS-62` guard,
fixed 2026-09-22, and the deliberate pre-rebuild `OPS-60` `skrf` red); all
four rungs solved and every imported gate green. Not re-run (2026-09-23
review). Successive changes by the review's hand arithmetic in the ledger
row: S₁₁ 2.815 → 1.497 % (1.88), S₂₁ 1.570 → 0.816 % (1.92), S₃₁ 1.362 →
0.831 % (1.64), with the finest cut-off rung's ≈ 10× C4 spread as the
caveat — the weekly applies the decision rule below. *(Was QUEUED
2026-09-18 for Tuesday 2026-09-22 02:00 —
`docs/testing/xl-queue.d/20260922-ANS-4-step3e.env`.)*

**Licence class:** priced-family variant of the `ANS-4-step2d` ledger row
(2026-09-10) — the same module, the same generator, **the same
`RUNGSPEC="0.015:1 0.015:2 0.0075:2 0.005:2"` as step 2d, frequency varied
only**: `FEM_EM_ANS4_FREQUENCY_HZ=64e6` (the knob `ANS-4-step3` /
`step3c` already ran green at `-n 16`). Every rung has a measured point at
this tier (step 2d: 116 085 / 116 085 / 281 728 / 592 744 cells, 0.14 /
0.74 / 1.80 / 3.79 M unknowns), so this is a variant, not a cost probe.

**Why (the §10 question it answers):** step 2d showed the 128 MHz degree-2
sequence *converging in h* (successive change 1.19 % → 0.62 %, ratio 1.90),
which is what lets the 128 MHz order-matched rung be read as a value rather
than a point. At 64 MHz only the two end rungs exist (`ANS-4-step3`:
`0.015:1` and `0.005:2`), so the 64 MHz order-matched figure the 09-19
weekly adjudicates has no *h*-convergence statement under it. This window
supplies the two missing degree-2 rungs at 64 MHz and the same successive-
change readout.

**Command:**

```
XL_CHUNK="ANS-4-step3e"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=64e6 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2 0.0075:2 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3e-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3e-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3e-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3e-raw.log; cat /workspace/logs/ans4-step3e-raw.log; exit \$rc'"
```

**Price, inherited from the family (scaling stated):** step 2d ran this
exact `RUNGSPEC` at `-n 16` in **7225 s, 290.2 GiB peak** *before* factor
reuse (four-drive times 12.0 / 169.0 / 1119.2 / 5653.5 s). With `PORT-19`
reuse in the module since 09-11 the finest rung's four drives measured
**1476.9 s** (`ANS-4-step3`), a factor 3.8; applying the same factor to the
two middle rungs gives ≈ 12 + 45 + 295 + 1480 s of solves plus four mesh
builds — **predicted 2 000–3 500 s**, memory set by the finest rung,
**270–295 GiB**. Timeout 14 400 s inside the 4 h window; 512 GiB limit.
Frequency does not change the factorisation's size.

**Readout (record, nothing asserted beyond the module's imported gates and
its own ladder controls):** public — the three C4 classes of `S` at 64 MHz
on the three degree-2 rungs, the successive changes and their ratio beside
step 2d's 128 MHz 1.19 % → 0.62 % (1.90), every imported `PORT-11` gate per
rung; private (`docs/private/`) — nothing new is needed: the AED column is
the one step 3 is read against. **Decision rule for the weekly:** successive
degree-2 changes falling with ratio ≳ 1.5 and the last one below the
64 MHz two-rung move's smallest class (step 3, 2.60 %) ⇒ step 3's `0.005:2`
figure is an *h*-converged value at the stated successive change and the
adjudication may quote it as such; a non-monotone or non-falling sequence
⇒ the 64 MHz order-matched figure is a point, not a value, and the row
says so. No band moves. The module's 128 MHz record assert is skipped at
`-n 16` by design (`OPS-41`). **Trap already paid for:** step 2d's single
failure was the ladder's refinement control comparing two rungs that share
one mesh (`0.015:1` vs `0.015:2`), fixed 2026-09-10; this spec has run
green in that module only at 128 MHz, so a red on a ladder-control test at
64 MHz is a module finding for a §9 item, not a physics reading.

### 8. `xl` — `ANS-4` step 3f: the 10 MHz degree-2 rung on the C4-congruent cut (daily licence: priced-family variant, 2026-09-18 review)

**Status:** RUN `20260923T070007Z_ANS-4-step3f.log` — 18 passed, Status 0,
1702 s, 280.28 GiB; degree-2 class spreads with the cut on **0.0303 / 0.0158
/ 0.0100 %** against 3c's 0.4451 / 0.2199 / 0.2186 % — the decision rule's
**collapse branch** by the public numbers (the weekly rules). Two-rung move
4.4499 / 1.2966 / 1.2187 %. Follow-up written by the 2026-09-23 review under
the daily licence: entry 12 (the 10 MHz *h*-ladder with the cut on). *(Was
QUEUED 2026-09-18 for Wednesday 2026-09-23 02:00 —
`docs/testing/xl-queue.d/20260923-ANS-4-step3f.env`.)*

**Licence class:** priced-family variant of `ANS-4-step3c` (2026-09-18,
1603 s / 284.8 GiB) and `ANS-4-step3b` (2026-09-17, 1539 s / 273.1 GiB) —
same module, same `RUNGSPEC`, same frequency as 3c (10 MHz); **one named
env knob varied**: `FEM_EM_ANS4_STEP2_C4_CONGRUENT=1`, the knob 3b already
ran green at this rung (592 550 cells, 3 827 268 dofs).

**Why (the status it can move):** step 3c's degree-2 rung passed every
imported `PORT-11` gate, but its **self-class C4 spread read 0.4451 %
against the 0.5 % band** (`20260918T070008Z_ANS-4-step3c.log:4624`) — 89 %
of the band, where the identical mesh at 64 MHz read 0.1985 % and the
congruent-cut rung at 128 MHz read 0.0056 %. One frequency pair cannot say
whether that is the cut-off mesh's asymmetry growing as f falls or a
low-frequency property of the degree-2 solve. This window changes the one
variable 3b changed.

**Command:**

```
XL_CHUNK="ANS-4-step3f"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=10e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3f-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3f-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3f-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3f-raw.log; cat /workspace/logs/ans4-step3f-raw.log; exit \$rc'"
```

**Price, inherited:** three measured windows on this `RUNGSPEC` at `-n 16`
— 1658 / 1539 / 1603 s, 280.8 / 273.1 / 284.8 GiB (all with a warm FFCx
cache; the named cache volume survives the pre-window recreate). Prediction
**1 500–1 800 s, 270–290 GiB**. Timeout 14 400 s; 512 GiB limit.

**Readout (record; the imported gates are the module's own and stay
unmoved):** the degree-2 rung's three class spreads at 10 MHz with the cut
on, beside 3c's 0.4451 / 0.2199 / 0.2186 % (cut off); the two-rung move on
the driven column beside 3c's 4.41 / 1.35 / 1.08 %. **Decision rule for the
weekly:** self spread collapsing by an order (as 128 MHz did, to the 1e-2 %
class) ⇒ the 0.4451 % is the cut-off port-sheet asymmetry, `GEO-32`'s row
gains the degree-2 10 MHz reading and the congruent cut is the rung to
quote at low frequency; self spread staying in the 0.4 % class ⇒ a
frequency trend in the degree-2 C4 identity that would breach the band
below 10 MHz — a known-issues entry and a `PORT` systematics item, **never
a band change**. A red `PORT-11` gate in this window (spread > 0.5 %) is
the same finding, stronger, and is reported as measured.

### 9. `xxl` — `WF-7` step 0d: the F-human degree-2 solve at the full 32-port set (daily licence: priced-family variant, 2026-09-19 review)

**Status:** RUN `20260926T070009Z_WF-7-step0d.log` — `[capture] rc=0`,
Status 0, 646 s, 103.06 GiB peak; full degree-2 32×32 reciprocity 9.715e-14
(asserted); terminal residual 1.08–1.18e-2 on every column. Readout and
decision-rule reading in `xxl-ledger.md` (2026-09-26 03:00 review). *(Was:
QUEUED 2026-09-19 for Saturday 2026-09-26 02:00 —
`docs/testing/xxl-queue.d/20260926-WF-7-step0d.env`.)*

**Licence class:** priced-family variant of the `WF-7-step0b` `xxl` ledger
row (2026-09-19) — the same probe, the same mesh (507 266 cells), the same
order (degree 2, measured at this tier: 3 256 418 unknowns, 347.56 s,
107.9 GiB summed, `memory.peak` 106.1 GiB), the same 64 MHz; **port set
varied only**: `FEM_EM_WF7_PORTS=all`, the knob `a267c5a` landed and entry 6
runs at degree 1 on 09-21. The two knobs are independent in the probe
(`wf7_step0_f_human_cost.py:117`, `:123`) but **have never run together** —
a red at the knob interaction is a probe finding for a §9 item, not a
physics reading. The weekly may replace or reorder this entry by renaming
or deleting the queue file; nothing depends on it.

**Why (the §10 question it answers):** step 0b priced *one* degree-2 column
at human scale and read a 5.50 % order move on it. If the weekly's branch
(a) holds (degree 2 is the production order at human scale), the quantity
§10 needs is the price of the **whole** degree-2 32×32 — one factorisation
plus 31 back-substitutions — and whether the human-scale degree-2 matrix
keeps `PORT-13`'s identities. It also answers, for free, whether step 0b's
degree-2 power-accounting residual (1.104e-02, 2.6× the degree-1 one) is a
property of port P17 or of every column.

**Command:**

```
XL_CHUNK="WF-7-step0d"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/wf7-step0d-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF7_DEGREE=2 FEM_EM_WF7_PORTS=all timeout -k 60 14400 mpiexec -n 16 python3 scripts/probes/wf7_step0_f_human_cost.py >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, inherited from the family (scaling stated):** step 0b's degree-2
leg — 133 s mesh, 347.56 s first drive (factorisation 242.5 s), MUMPS
forward/backward solve driver **1.846 s** on the held factor
(`20260919T070008Z_WF-7-step0b.log:21167`). Thirty-one further columns at a
few seconds each (the degree-1 two-drive window read 0.59 s per held
column against a 0.20 s solve driver, a factor ≈ 3 of overhead) ⇒
**predicted 10–20 min, 105–130 GiB** (the factor is held once; column fields
are dropped after each `S` column is read). Timeout 14 400 s inside the 8 h
window; 754 GiB limit. The probe's own degree-2 brackets (150–350 GiB,
10–60 min *per drive*) are step 0b's stale prediction and, unless `OPS-52`
has re-dated them by then, will print OUTSIDE — a label, not a finding
(known-issues 2026-09-19).

**Readout (record — nothing asserted beyond the imported cell band):**
unknowns, first-drive and per-held-column solve times, cumulative wall
clock, `ru_maxrss` per rank / summed, `memory.peak`; the degree-2 32×32's
`_reciprocity_ratio`, `σ_max` and class spreads **printed** beside
`RECIPROCITY_BAND` / `COLUMN_PASSIVITY_CEILING` / `OPPOSITE_SPREAD_BAND`,
next to entry 6's degree-1 values; the per-column power-accounting
residual. **Decision rule for the weekly:** a full degree-2 set inside one
window at this price ⇒ the human-scale degree-2 S sweep is an ordinary
`xl`-sized job (106 GiB fits the 512 GiB service) and §10 dates it as such;
a reciprocity or passivity miss at degree 2 that entry 6 does not show at
degree 1 ⇒ a finding for a `PORT` chunk; residuals near 1e-2 on every
column ⇒ the accounting's degree-2 quadrature or the sheet term is the
suspect, a `POST` item — **never a band change**.

### 10. `xl` — `ANS-2` step 4: the phantom h-halving, four drives (operator session 2026-09-19; daily-licence class: cost probe, price measured by the step's own probe)

**Status:** RUN `20260925T070008Z_ANS-2-step4.log` — Status 0, 1315 s,
39.91 GiB, orphans 0 / 0: C4 driven-point spread 15.730 → 4.528 %, phantom
power moves ≤ 0.30 % ⇒ **branch (a)**, ruled by the 2026-09-25 review (ledger
row; §7 `ANS-2`). `metrics_h0.00125.json` is **not committed**, despite
"Outputs to commit" below. The `OPS-59` artifact pin
(`tests/unit/test_doc_reference_exit_codes.py:463`) lists exactly one
`metrics.json` per `ans:` case, so a second tracked JSON would turn it red.
The review moved the file beside the XDMF under the gitignored
`paraview_output/`; the log carries every readout figure. *(Was QUEUED
2026-09-20 for Friday 2026-09-25 02:00 —
`docs/testing/xl-queue.d/20260925-ANS-2-step4.env`; before that READY.)* `ANS-2` step 4a landed 2026-09-19 (§9 item 9): the example
reads `FEM_EM_ANS2_PHANTOM_RESOLUTION`, tags `metrics_h<res>.json` and the
XDMF, skips the two 0.0025-rung mesh-record assertions **with a printed
`[ANS-2 step 4] MESH-RECORD ASSERTIONS SKIPPED` line** while every other
imported band stays asserted, leaves `COMPARISON.md` /
`COMPARISON_private.md` alone, and prints the `[ANS-2 step 4]` readout
(driven-point SAR, C4 spread, phantom powers, each beside the tracked 0.0025
rung's own `metrics.json`). Anchors measured green — unset vs the tracked
`metrics.json` and knob-`0.0025` vs unset both digit-identical at
`EXACT_IDENTITY_RTOL` on every primary leaf (worst 8.19e-11 / 8.17e-11),
negative control at 0.004 m prints the skip line and exits 0
(`20260919T141006Z_ANS-2-step4a-unset-control.log`,
`20260919T142050Z_ANS-2-step4a-knob-identity.log`). *The next review only
queues this entry — the first free `xl` night (Thu 2026-09-24 02:00 as of
writing).*

**Why `xl`:** the step's pre-registered cost probe fired its own STOP rule
(`32f4eae`, `20260919T093241Z_ANS-2.log`): at `phantom_resolution` =
0.00125 m the fixture meshes to **719 769 cells (452 228 phantom) / 845 188
unknowns**, mesh **359.8 s**, one drive **188.4 s**, `ru_maxrss` **40.89 GiB
summed** at `-n 8` — past the item's 15 min / 40 GiB stop and past §5.1's
20-min cap for four drives. §5.1's rule is met: the price is *measured*, on
this mesh, at this width.

**What runs:** the `ANS-2` example at the halved rung, all four single-port
drives, `-n 8` (the measured width; 16 ranks is unmeasured on this mesh and
not needed), durable capture per §5.1. **Outputs to commit:**
`metrics_h0.00125.json` and
`paraview_output/ans2_birdcage_coil_driven_sar_10mhz_combined_h0p00125.xdmf`
(the XDMF tag writes the decimal point as `p` — `write_xdmf_with_tags` strips
everything after the last dot as a suffix, which silently collapsed every
rung onto one filename until step 4a measured it); `metrics.json` and
`COMPARISON.md` are untouched by an overridden run.

**Command (final once 4a lands; the knob name is the contract):**

```
XL_CHUNK="ANS-2-step4"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans2-step4-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS2_PHANTOM_RESOLUTION=0.00125 timeout -k 60 7200 mpiexec -n 8 python3 examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/02_birdcage_coil_driven_sar_10MHz.py >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, measured (the probe):** mesh 360 s once + 4 × 188 s ≈ **1 115 s** —
step 4a settles the open branch: **the example meshes once** for all four
drives (0.0025 rung, `-n 4`: 184.8 s wall of which 74.4 s is the four solves,
so a single mesh build precedes them; `20260919T142050Z_ANS-2-step4a-knob-identity.log`),
so the ≈ 2 200 s re-meshing arm is dead; ≈ 41 GiB summed at `-n 8` (the probe's figure; four drives hold one
field at a time). Timeout 7 200 s inside the 4 h window; 512 GiB limit.

**Readout (record — the imported bands asserted, the mesh-record bands
skipped with a printed line per 4a):** the four driven-point SAR values and
their C4 spread beside the 0.0025 rung's (15.7 %); the four phantom powers
beside the 0.0025 rung's and their relative change; the mass-averaged C4
pairs at the imported 5 % band; the halved rung's cell counts as a
candidate record. **Decision rule (pre-registered in
`docs/private/ans2-adjudication-2026-09-18.md` and §7 `ANS-2` step 4):**
(a) spread < 5 % and phantom power moves < 2 % ⇒ the residual against
HFSS is the feed; the verdict's pointwise band tightens to the measured
spread. (b) phantom power moves ≥ 2 % toward AED ⇒ mesh convergence is part
of the residual; a third rung is priced by a review. (c) spread does not
fall ⇒ not `h` — a point-evaluation defect, known-issues + `POST` chunk.
No band moves in the window; the review that finds the ledger row rules.

### 11. `xl` — `ANS-4` step 3g: the 10 MHz degree-2 *h*-ladder (weekly 2026-09-19, §10 chain F step F2)

**Status:** RUN `20260924T070008Z_ANS-4-step3g.log` — 18 passed, Status 0,
2084 s, 281.95 GiB (ledger row, 2026-09-25 review): gates green on every
rung; successive degree-2 changes 6.969 → 3.606 % / 1.813 → 0.934 % /
1.843 → 1.100 % (ratios 1.93 / 1.94 / 1.67) — falling, but the last change
is above 1.08 % on S₁₁ and S₃₁, so neither branch of the rule below is met
as written; the weekly rules, with 3h (entry 12) as the cut control. *(Was
QUEUED 2026-09-20 for Thursday 2026-09-24 02:00 —
`docs/testing/xl-queue.d/20260924-ANS-4-step3g.env` (ahead of entry 10 so
that a budget-denied night still leaves 3g run before the 09-26 weekly,
which reads it with 3e). Before that READY.)* No prerequisite; the knob, the module and the `RUNGSPEC`
have all run green at `-n 16` (`ANS-4-step3c` 2026-09-18 at 10 MHz, step 2d
2026-09-10 with this `RUNGSPEC`). The next daily review queues it for the
first free `xl` night.

**Why (the §10 question it answers):** the 2026-09-19 weekly ruled that 3c's
public rule fired — the 10 MHz degree 1 → 2 move (4.41 / 1.35 / 1.08 %) is in
the same class as the Larmor moves — and that the order-matched comparison
disagrees on the self class at 10 MHz (known-issues 2026-09-19; numbers
private). Both readings rest on **one** degree-2 rung at 10 MHz. Entry 7
(3e) supplies the missing degree-2 rungs at 64 MHz; this is the same ladder
at 10 MHz, where the disagreement is largest.

**Command:** entry 7's, with `FEM_EM_ANS4_FREQUENCY_HZ=10e6` and the raw-log
name changed — nothing else:

```
XL_CHUNK="ANS-4-step3g"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=10e6 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2 0.0075:2 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3g-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3g-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3g-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3g-raw.log; cat /workspace/logs/ans4-step3g-raw.log; exit \$rc'"
```

**Price, measured in the family:** the finest rung at 10 MHz is 3c's own
window — four drives 1419.4 s, `memory.peak` 284.8 GiB, 1603 s total; the two
middle rungs scale from step 2d by the reuse factor entry 7 states
(≈ 45 s + ≈ 295 s) plus four mesh builds ⇒ **predicted 2 000–3 500 s,
270–295 GiB**. Frequency does not change the factorisation's size. Timeout
14 400 s inside the 4 h window; 512 GiB limit.

**Readout (record; nothing asserted beyond the module's imported gates and
ladder controls):** public — the three C4 classes at 10 MHz on the three
degree-2 rungs, the successive changes and their ratio beside step 2d's
128 MHz 1.19 % → 0.62 % (1.90) and 3e's 64 MHz pair; every imported gate per
rung (watch the self-class spread: 0.4451 % of a 0.5 % band on 3c's rung —
a red there is a `GEO-31`-class mesh reading, and entry 8 (3f, the congruent
cut at 10 MHz) is its control, not a reason to touch the band). Private —
nothing new: the AED column is the one 3c was read against.

**Decision rule (pre-registered, PROJECT_PLAN §10 2026-09-19):** successive
degree-2 changes falling with ratio ≳ 1.5 **and** the last change below the
smallest class move (1.08 %) ⇒ 3c's `0.005:2` figure is an *h*-converged
value, the 10 MHz self-class disagreement is a model-level difference, and
with the same outcome on 3e the production default flips to degree 2 by
`TH-19` step 4 and `PORT-21` step 2 reads its table against a converged
residual. Non-monotone or non-falling ⇒ the 10 MHz order-matched figure is a
point, not a value; the default stays at degree 1; the ladder is the finding.
No band moves.

### 12. `xl` — `ANS-4` step 3h: the 10 MHz degree-2 *h*-ladder on the C4-congruent cut (daily licence: priced-family variant, 2026-09-23 review)

**Status:** RUN `20260927T070008Z_ANS-4-step3h.log` — 1994 s (ledger row;
reading by the 2026-09-27 03:00 review). *Was:* QUEUED 2026-09-23 for Sunday
2026-09-27 02:00 — `docs/testing/xl-queue.d/20260927-ANS-4-step3h.env`.

**Licence class:** priced-family variant of entry 11 (`ANS-4` step 3g, the
weekly's chain step F2) — the same module, the same
`RUNGSPEC="0.015:1 0.015:2 0.0075:2 0.005:2"`, the same 10 MHz; **one named
env knob varied**: `FEM_EM_ANS4_STEP2_C4_CONGRUENT=1`, the knob 3b / 3d / 3f
ran green at `-n 16` on this module (3f at this frequency, 2026-09-23).

**Why (the §10 question it answers):** F2 asks whether the 10 MHz degree-2
rung is *h*-converged. Entry 8 (3f) has since come back in its **collapse
branch** — with the cut on, the degree-2 self-class C4 spread at 10 MHz falls
from 3c's 0.4451 % to 0.0303 % — so the cut-off mesh entry 11 runs carries a
≈ 0.4 % port-sheet asymmetry on its finest rung, the same size as the
successive changes F2 is trying to read (3e's 64 MHz cut-off ladder shows the
mechanism: its finest rung's spread is ≈ 10× the middle rungs', which leaves
S₃₁'s ratio anywhere in ≈ 1.3–2.2). Entry 8's own rule names the congruent
cut as the rung to quote at low frequency; this is F2's ladder on that mesh.
Entry 11 is not replaced — the pair separates *h* from the cut.

**Command** (entry 11's, with the knob added exactly where entry 8 carries it
and the raw-log name changed — nothing else):

```
XL_CHUNK="ANS-4-step3h"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=10e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2 0.0075:2 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3h-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3h-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3h-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3h-raw.log; cat /workspace/logs/ans4-step3h-raw.log; exit \$rc'"
```

**Price, measured in the family (scaling stated):** this `RUNGSPEC` at
`-n 16` ran in **2206 s, 277.44 GiB** (3e, cut off, 64 MHz; four-drive
times 4.3 / 52.0 / 309.6 / 1502.0 s); the finest rung **with the cut on at
10 MHz** is 3f's own — 592 550 cells, 3 827 268 dofs, four drives 1476.9 s,
280.28 GiB. The cut-on 0.0075 rung's cell count is unmeasured (cut off:
281 728); memory is set by the finest rung. Frequency does not change the
factorisation's size ⇒ **predicted 2 100–2 700 s, 275–290 GiB**. Timeout
14 400 s inside the 4 h window; 512 GiB limit.

**Readout (record; the module's imported gates and ladder controls, nothing
else asserted):** the three C4 classes and the imported `PORT-11` gates on
each of the three degree-2 rungs with the cut on; the printed driven-column
entries per rung. **The module does not print successive changes** — the
review that finds the ledger row computes them by hand (module's move
convention, coarser rung as reference) as the 2026-09-23 review did for 3e,
beside entry 11's cut-off ladder. The finest rung reproduces 3f's (cells
592 550, spreads 0.0303 / 0.0158 / 0.0100 %) — printed, not asserted.
Private: nothing new.

**Decision rule for the weekly:** successive degree-2 changes falling with
ratio ≳ 1.5 on **all three** classes **and** the last change below 3f's
smallest two-rung class move (1.2187 %) ⇒ 3f's cut-on `0.005:2` figure is an
*h*-converged value at 10 MHz and F2's 10 MHz answer is read on this ladder;
entry 11 agreeing within the cut-off rung's own spreads ⇒ the cut does not
move the *h*-statement; disagreeing ⇒ entry 11's reading is the cut's, and
this entry's is the one F2 records. Non-monotone or non-falling ⇒ the 10 MHz
order-matched figure is a point, not a value. No band moves. A red on a
ladder-control test is a module finding for a §9 item, not physics.

### 13. `xl` — `ANS-4` step 3i: the 64 MHz degree-2 *h*-ladder on the C4-congruent cut (daily licence: priced-family variant, 2026-09-23 review)

**Status:** RUN `20260928T070008Z_ANS-4-step3i.log` (Status 0, 2241 s,
282.44 GiB; ledger row read by the 2026-09-28 03:00 review). *Was:* QUEUED
2026-09-23 for Monday 2026-09-28 02:00 —
`docs/testing/xl-queue.d/20260928-ANS-4-step3i.env`.

**Licence class:** priced-family variant of entry 7 (`ANS-4` step 3e,
**RUN 2026-09-22, 2206 s / 277.44 GiB**) — the same module, `RUNGSPEC` and
64 MHz; **one named env knob varied**: `FEM_EM_ANS4_STEP2_C4_CONGRUENT=1`,
run green on this module's finest rung at 64 MHz by 3d (2026-09-20).

**Why (the §10 question it answers):** 3e's successive changes (ledger row,
review arithmetic) are S₁₁ 2.815 → 1.497 %, S₂₁ 1.570 → 0.816 %, S₃₁ 1.362 →
0.831 %, ratios 1.88 / 1.92 / 1.64 — but its finest rung is the cut-off mesh
whose C4 spread (0.17–0.20 %) is ≈ 10× the middle rungs', so a last change of
≈ 0.8 % read on one C4 member carries ≈ ±0.2 pp of member choice and S₃₁'s
ratio straddles the 1.5 threshold. With the cut on the 64 MHz finest rung
reads 0.0132 / 0.0126 / 0.0081 % (3d). This window removes the ambiguity
from the F2 answer the weekly rules on at 64 MHz.

**Command** (entry 7's, with the knob added where entry 8 carries it and the
raw-log name changed — nothing else):

```
XL_CHUNK="ANS-4-step3i"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=64e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2 0.0075:2 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3i-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3i-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3i-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3i-raw.log; cat /workspace/logs/ans4-step3i-raw.log; exit \$rc'"
```

**Price, measured in the family:** 3e — this `RUNGSPEC`, cut off, 64 MHz,
`-n 16` — **2206 s, 277.44 GiB**; 3d — the cut-on finest rung at 64 MHz —
1619 s, 281.1 GiB (four drives 1425.5 s) ⇒ **predicted 2 100–2 700 s,
275–290 GiB**. Timeout 14 400 s; 512 GiB limit.

**Readout (record):** as entry 12, at 64 MHz, beside 3e's cut-off ladder;
the finest rung reproduces 3d's (592 550 cells, 0.0132 / 0.0126 / 0.0081 %)
— printed, not asserted. Successive changes by the review's hand
arithmetic, as for 3e. Private: nothing new.

**Decision rule for the weekly:** entry 7's, unchanged — ratio ≳ 1.5 on all
three classes and the last change below step 3's smallest class move
(2.60 %) ⇒ the 64 MHz order-matched figure is an *h*-converged value; now
read on the mesh whose finest rung does not carry the 10× C4 spread. 3e and
this entry agreeing within 3e's finest-rung spread ⇒ 3e's answer stands;
disagreeing ⇒ this entry's is the one F2 records. No band moves.

### 14. `xl` — `ANS-6` step 2c: the hole-mesh degree-2 rung at h = 0.005, 64 MHz (daily licence: cost probe, 2026-09-23 review; chain T9b's second half)

**Status:** RUN `20260929T070008Z_ANS-6-step2c.log` — `[capture] rc=0`,
Status 0, 4110 s, every imported gate green (`:1037–1041`); launcher
orphans 0 / 0 (`:1054–1059`). **Readout against the brackets (2026-09-29
10:30 review): INSIDE on every one** — 523 816 cells (350 k–600 k),
3 466 952 unknowns (2.2–3.9 M), factorisations 1 461.2 / 1 302.3 /
1 236.9 s (600–1 500 s), total 4 106.1 s (2 000–5 000 s), `memory.peak`
266.79 GiB (286 461 792 256 B; 170–300 GiB), summed `ru_maxrss` 295.22 GiB
(`:1040`, `:1044–1046`). Decision rule: the matched hole rung is an ordinary
`xl` window ⇒ entries 17 and 18 marked READY and queued. Printed, not
gated: the step's predicted σ-vs-PEC bound (≤ 1e-4) is **not met** at
4.008e-04 max entry (`:1033`); the class moves against degree 1 read
4.48 / 2.00 / 2.95 % (Cu), predicted 3–7 % (`:1038`). *Was:* QUEUED
2026-09-25 for Tuesday 2026-09-29 02:00 —
`docs/testing/xl-queue.d/20260929-ANS-6-step2c.env` (command line
byte-identical to the block below). *Before that* READY (2026-09-23, 04:30 CDT slot) — §9 item 35 (`ANS-6` step 2b,
the `FEM_EM_ANS6_RESOLUTION` knob) landed: control at h = 0.015 reproduces
the tracked 64 MHz S leaves to 2.647e-15 on 80 181 cells / 111 121 unknowns
(`20260923T094235Z_ANS-6-step2b-a-h0p015-d1-64MHz.log:940–944`); the knob
reaches the mesh at h = 0.0075 (226 642 cells / 294 763 unknowns, 147.4 s,
6.39 GiB summed RSS at degree 1, `-n 4`,
`20260923T093956Z_ANS-6-step2b-b-h0p0075-d1-64MHz.log:950–953`). The next
review queues it on the first free `xl` night (Tue 2026-09-29 as of writing
— 09-24 … 09-28 are taken). *Was:* PENDING PREREQUISITE on item 35.

**Licence class:** cost probe — the first window on an unpriced case (no
point exists on the *hole* mesh finer than h = 0.015). Chain T9b (§10,
2026-09-19): "an `xl` cost probe on the hole mesh at h = 0.005 degree 2".
The readout *is* the price; fixed brackets, no gate, no band.

**Why:** the `ANS-6` preliminary reading (§7 row, 2026-09-20) attributes
`ANS-4`'s self-class offset to the lumped-port model / our degree-1 mesh and
clears the resolved-conductor interior; the hole mesh removes the interior
from both codes (chain F5). The order- and *h*-matched hole rung is what the
09-26 weekly's adjudication of `ANS-6` needs next, and nothing prices it.

**Command** (the knob name is the contract; final when item 35 lands):

```
XL_CHUNK="ANS-6-step2c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans6-step2c-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS6_DEGREE=2 FEM_EM_ANS6_FREQ_MHZ=64 FEM_EM_ANS6_RESOLUTION=0.005 timeout -k 60 14400 mpiexec -n 16 python3 examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz/06_copper_birdcage_four_port_10_64_128MHz.py >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, predicted (brackets for the readout):** the knob route runs three
factorisations (PEC sweep, Cu sweep, P1). At h = 0.015 degree 2 the hole
mesh (80 181 cells) took 61.9 / 53.9 / 44.9 s, 183.9 s, 13.78 GiB summed at
`-n 8` (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:924–942`). The
nearest measured h = 0.005 degree-2 point is `ANS-4`'s *solid-coil* mesh:
592 744 cells, 3 828 494 unknowns, factorisation 1 332 s, 277 GiB at `-n 16`
(3e). The hole mesh drops the conductor interior (80 181 / 116 085 = 0.69 at
h = 0.015) ⇒ **brackets: 350 k–600 k cells, 2.2–3.9 M unknowns, 600–1 500 s
per factorisation, total 2 000–5 000 s, 170–300 GiB `memory.peak`**, printed
INSIDE / OUTSIDE by the review, never asserted. Timeout 14 400 s inside the
4 h window; 512 GiB limit.

**Readout (record — the script's own imported gates asserted as in step 2,
no band):** cells, unknowns (ZMUMPS `N`), mesh and per-factorisation times,
`memory.peak`, summed `ru_maxrss`; both columns' reciprocity / σ_max / C4
spreads and the copper surface-loss identity (the step-2 gates); the class
moves against the h = 0.015 degree-2 point, printed. The untracked
`metrics_degree2_…` file is not committed. **Decision rule:** inside the
brackets ⇒ the matched hole rung is an ordinary `xl` window and the weekly
may commission it at 10 and 128 MHz; outside above ⇒ the price is the
finding; a gate red ⇒ a finding for a `TH-15` / `TH-14` §9 item, never a
band change.

### 15. `xl` — `ANS-4` step 3j: the 128 MHz degree-2 *h*-ladder on the C4-congruent cut (daily licence: priced-family variant, 2026-09-25 review)

**Status:** RUN `20260930T070008Z_ANS-4-step3j.log` — 18 passed, 2 skipped
(the expected frequency-knob skip and the conductor-factor test, env unset),
`[capture] rc=0`, Status 0, 2315 s, 278.59 GiB; successive degree-2 changes
S₁₁ 1.200 → 0.626 % (ratio 1.92), finest rung = 3b's (ledger row, filled by
the 2026-09-30 03:00 review). The reading against the decision rule is the
weekly's. *(Was: QUEUED 2026-09-25 for Wednesday 2026-09-30 02:00 —
`docs/testing/xl-queue.d/20260930-ANS-4-step3j.env`.)*

**Licence class:** priced-family variant of entry 13 (`ANS-4` step 3i) —
the same module, `RUNGSPEC="0.015:1 0.015:2 0.0075:2 0.005:2"` and
congruent cut (`FEM_EM_ANS4_STEP2_C4_CONGRUENT=1`); **one knob varied: the
frequency**, by leaving `FEM_EM_ANS4_FREQUENCY_HZ` unset (the module's
default, 128 MHz — how entry 3 (3b) and step 2d ran it). Both halves are
measured at `-n 16`: this `RUNGSPEC` at 128 MHz cut off (step 2d, 2026-09-10)
and the cut-on finest rung at 128 MHz (3b, 2026-09-17). Written because the
`xl` floor was one short after entry 14 was queued (§6b rule 4); the weekly's
list of preferred variants names "the `ANS-4` degree-2 rung at the other
Larmor frequency".

**Why (the §10 question it answers):** the 09-26 weekly's production-order
ruling (chain T6, `TH-19` step 4) re-records every degree-1 S record beside
its degree-2 value, at 10, 64 and 128 MHz. With 3h (10 MHz) and 3i (64 MHz)
queued, 128 MHz is the only frequency whose degree-2 *h*-ladder rests on the
**cut-off** mesh (step 2d: successive change 1.19 → 0.62 %, ratio 1.90). Its
finest rung carried the sheet-triangulation spread that 3b later removed
(0.0056 / 0.0088 / 0.0070 % with the cut on). 3g (entry 11, 2026-09-24) has
just shown that a ladder can fall with ratio ≈ 1.9 while its last change is
still ≈ 3.6 % on the self class. This window puts all three frequencies on one
mesh family, so T6 reads its ladders on the same cut.

**Command** (entry 13's with `FEM_EM_ANS4_FREQUENCY_HZ=64e6` removed and the
raw-log name changed — nothing else):

```
XL_CHUNK="ANS-4-step3j"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2 0.0075:2 0.005:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step3j-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step3j-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step3j-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step3j-raw.log; cat /workspace/logs/ans4-step3j-raw.log; exit \$rc'"
```

**Price, measured in the family:** this `RUNGSPEC` cut off at `-n 16` —
2206 s / 277.44 GiB (3e, 64 MHz) and 2084 s / 281.95 GiB (3g, 10 MHz), both
with factor reuse. Step 2d's 7225 s at 128 MHz predates held
back-substitution and is not the price. The cut-on finest rung at 128 MHz is
3b's own, four drives 1361.4 s, 273.1 GiB. Frequency does not change the
factorisation's size ⇒ **predicted 2 000–2 700 s, 270–290 GiB**. Timeout
14 400 s inside the 4 h window; 512 GiB limit.

**Readout (record; the module's imported gates and ladder controls, nothing
else asserted):** as entries 12 and 13, at 128 MHz, beside step 2d's cut-off
ladder. The finest rung reproduces 3b's (592 550 cells, 0.0056 / 0.0088 /
0.0070 %), printed, not asserted. **Expected skip:** the frequency-knob test
takes its `is_default` branch and skips at `-n 16` against the `-n 8`
record width (`OPS-41`), exactly as 3b did (15 passed / 1 skipped). A skip
there is not a red. Successive changes are the review's hand arithmetic, as
for 3e and 3g. Private: nothing new; the 128 MHz AGREE was read on 3b.

**Decision rule for the weekly:** step 2d's reading (ratio ≳ 1.5, the last
change below step 2b's smallest 128 MHz degree 1 → 2 class move, 5.3841 %,
ledger 2026-09-09 — the analogue of entry 7's 2.60 % and entry 11's 1.08 %)
either reproduces on the cut —
the 128 MHz degree-2 figure is an *h*-converged value on both meshes — or it
does not, in which case this entry's is the one T6 records. No band moves. A
red on a ladder-control test is a module finding for a §9 item, not physics.

### 16. `xl` — `ANS-6` step 2d: the hole-mesh degree-2 rung at h = 0.0075, 64 MHz (weekly 2026-09-26, §10 chain T9b)

**Status:** RUN `20261001T070008Z_ANS-6-step2d.log` — `[capture] rc=0`,
Status 0, 818 s, `memory.peak` 79 991 607 296 B (74.5 GiB), orphans none
(log footer; marked by the 2026-10-01 03:00 review). Readout and the
pre-registered row (a) rule against the private column: the weekly's.
*Was:* QUEUED 2026-09-27 for Thursday 2026-10-01 02:00 —
`docs/testing/xl-queue.d/20261001-ANS-6-step2d.env` (command line
byte-identical to the block below; the first free `xl` night, 5 of 6 charged
at its cutoff). *Was:* READY (2026-09-26) — every knob landed with `ANS-6`
step 2b (§9 item 35, 2026-09-23). Independent of entry 14; may run before or
after it.

**Why:** `ANS-6`'s absolute S classes were ruled INCONCLUSIVE at degree 1
because the gap to AED is below our own degree 1 → 2 move on the coarse hole
mesh. Entry 14 prices the finest rung; this is the middle one, so the
hole-mesh degree-2 sequence at 64 MHz has three rungs (h = 0.015 on record,
0.0075 here, 0.005 from entry 14) and can be read under §10's 2026-09-26
criterion.

**Command** (entry 14's with the resolution and the raw-log name changed —
nothing else):

```
XL_CHUNK="ANS-6-step2d"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans6-step2d-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS6_DEGREE=2 FEM_EM_ANS6_FREQ_MHZ=64 FEM_EM_ANS6_RESOLUTION=0.0075 timeout -k 60 14400 mpiexec -n 16 python3 examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz/06_copper_birdcage_four_port_10_64_128MHz.py >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, measured in the family (scaling stated):** degree 2 at h = 0.015 on
this mesh — 80 181 cells, three factorisations, 183.9 s, 13.78 GiB summed at
`-n 8` (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:924–942`); the
h = 0.0075 mesh is measured at 226 642 cells (2.83×,
`20260923T093956Z_…-b-h0p0075-d1-64MHz.log:950–953`). Memory ∝ N^1.33 and
time ∝ N² from `ANS-4`'s four priced degree-2 rungs ⇒ **predicted 40–80 GiB,
1 200–2 400 s**. Timeout 14 400 s; 512 GiB limit.

**Readout (record — the script's imported gates asserted, no band):** cells,
unknowns, times, `memory.peak`; both columns' reciprocity / σ_max / C4
spreads and the copper identity; the three S classes on both columns,
printed. Private: the review reads the classes against the order-matched
AED column.

**Decision rule:** with entry 14, the last successive degree-2 change ≥ 3×
below the gap to AED on a class ⇒ that class is read as converged and
`ANS-6` row (a) is ruled AGREE or DISAGREE at the stated level; otherwise
row (a) stays INCONCLUSIVE and the ladder is the finding. A gate red is a
`TH-14` / `TH-15` §9 item, never a band change.

### 17. `xl` — `ANS-6` step 2e: the hole-mesh degree-2 rung at h = 0.005, 10 MHz (weekly 2026-09-26, §10 chain T9b)

**Status:** SUPERSEDED 2026-09-30 (weekly) — not run; the queue file was
removed in the weekly's commit. Reason: an *uncorrected* fine rung at
10 MHz repeats a known answer (entry 14 read 3.1 / 1.8 / 1.5 % at 64 MHz,
the same on Cu and PEC; 3h / 3i / 3j each reproduced a cut-off twin), and
the night goes to the κ-corrected rung (entry 22). **Revived as written
only if entries 22 and 23 reject the κ reading** (§10 rulings 2026-09-30).
*Was:* QUEUED 2026-09-29 (10:30 review) for Monday 2026-10-05 02:00 —
`docs/testing/xl-queue.d/20261005-ANS-6-step2e.env`. *Before that* READY
(2026-09-29 10:30 review): entry 14 landed inside every bracket with
`[capture] rc=0` (`20260929T070008Z_ANS-6-step2c.log:1037–1046`). *Was:*
PENDING PREREQUISITE — entry 14 (`ANS-6` step 2c, 2026-09-29)
landing **inside its brackets** with `[capture] rc=0`. The daily review marks
this READY on that condition and on no other; outside the brackets it is
superseded and the weekly re-prices.

**Command:** entry 14's, with `FEM_EM_ANS6_FREQ_MHZ=10`, `XL_CHUNK="ANS-6-step2e"`
and the raw-log name `ans6-step2e-raw.log` — nothing else changed.

**Price:** entry 14's measured one (frequency does not change the
factorisation's size). **Readout and decision rule:** as entry 16, at
10 MHz, where `ANS-4`'s order-matched self-class gap is largest — the PEC
column at this rung is the cleanest reading of the port model available.

### 18. `xl` — `ANS-6` step 2f: the hole-mesh degree-2 rung at h = 0.005, 128 MHz (weekly 2026-09-26, §10 chain T9b)

**Status:** SUPERSEDED 2026-09-30 (weekly) — not run; queue file removed
in the weekly's commit; reason and revival condition as entry 17. *Was:*
QUEUED 2026-09-29 (10:30 review) for Tuesday 2026-10-06 02:00 —
`docs/testing/xl-queue.d/20261006-ANS-6-step2f.env`. *Before that* READY
on entry 17's condition (met). *Was:* PENDING PREREQUISITE — as entry 17.

**Command:** entry 14's, with `FEM_EM_ANS6_FREQ_MHZ=128`,
`XL_CHUNK="ANS-6-step2f"` and the raw-log name `ans6-step2f-raw.log`.

**Price, readout, decision rule:** as entry 17, at 128 MHz — the control:
`ANS-4` agrees there at 1.1 % order-matched, so a hole-mesh gap much larger
than that at 128 MHz would be a finding about the hole route, not the port.

### 19. `xl` — `PORT-21` step 3b: the κ-corrected sheet width at degree 2 on the finest congruent rung (weekly 2026-09-26, §10 chain F7)

**Status:** SUPERSEDED 2026-09-30 (weekly) by **entry 23** — the same
command re-issued for the `xxl` night of 2026-10-03 at its measured price
(every `xl` night to 10-06 is taken and the 10-03 weekly needs the
reading), with the decision rule re-stated: at degree 2 on the ×1 cut (a)'s
move is 5.4× its degree-1 value, a case none of the three branches below
covers. *Was:* PENDING PREREQUISITE — `PORT-21` step 3 (chain F6, an ordinary
§9 item): the degree and resolution knobs on
`tests/validation/test_port21_feed_sensitivity.py`. **The knob names below
are the contract**; if the item lands other names this entry is re-issued,
not edited.

**Measured price beside it (`PORT-21` step 3, 2026-09-27, left PENDING
for the weekly):** F6 landed the three knob names exactly as the command
below spells them (space-separated frequencies accepted). W3 (degree 2,
10 MHz, `control_prime a`, ×1 cut, 746 144 dofs, `-n 8`,
`20260927T110931Z_PORT-21-step3-W3.log`) measured **three factorisations per
frequency, not two**: control′'s four-drive sweep 72.6 s, the in-run κ
solve ≈ 65.8 s (its own factorisation), (a)'s four-drive sweep 69.7 s, and
control′ re-meshes at every frequency (24.8 s here). Memory does **not**
stack: live RSS returns to ≈ 3 GiB after every sweep, peak sum
18.3 → 19.2 GiB across phases, so one factor at a time (the 280 GiB family
figure stands). Scaled on the 1 476.9 s family sweep (κ solve ≈ 0.9 of a
sweep) the "10 64" command is ≈ 5.8 sweep-equivalents + 2 mesh builds ≈
**8 600–9 500 s**, outside the 5 500–7 000 s bracket above and inside the
14 400 s timeout. Also printed there: degree-2 κ = 5.209e-02 vs 1.059e-02
at degree 1, (a)'s self move 1.846e-01 vs 3.411e-02.

**Why:** variant (a) has the frequency shape of the order-matched
self-class gap and points toward AED, but it was measured at degree 1 on the
×1 mesh. Whether its move survives at degree 2 on the finest rung is the
one unmeasured step in the §10 2026-09-26 reading.

**Command** (final when F6 lands):

```
XL_CHUNK="PORT-21-step3b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/port21-step3b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_PORT21_FREQUENCIES_MHZ=\"10 64\" FEM_EM_PORT21_DEGREE=2 FEM_EM_PORT21_RESOLUTION=0.005 FEM_EM_PORT21_VARIANTS=\"control_prime a\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_port21_feed_sensitivity.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, measured in the family:** the same mesh and order as `ANS-4` step
3f — 592 550 cells, 3 827 268 unknowns, one factorisation + four drives
1 476.9 s, 280.28 GiB at `-n 16`. Two variants × two frequencies = four
factorisations on one mesh build ⇒ **predicted 5 500–7 000 s, 275–295 GiB**.
Timeout 14 400 s; 512 GiB limit. If F6's heavy-tier price at degree 2 on
the ×1 cut contradicts the family scaling, the daily does not mark this
READY and returns it.

**Readout (record; the module's imported gates asserted on every run):**
the self / adjacent / opposite move of (a) against control′ at degree 2,
beside the degree-1 3.41 / 1.55 % (10 / 64 MHz). Private: both variants
against the order-matched AED column.

**Decision rule:** (a)'s self-class move at degree 2 within a factor 1.5 of
its degree-1 value **and** toward AED ⇒ the sheet-width convention is a
named, measured part of the gap and `PORT-14`'s κ correction becomes the
convention quoted against AED (a `PORT-21` step 4 re-record, no band
moved); the move collapsing at degree 2 ⇒ (a) was a coarse-mesh effect and
the gap is truncation alone (F8 decides); the move away from AED ⇒ the
hypothesis is dead and is written as one.

### 20. `xl` — `ANS-4` step 4b: the first degree-2 rung that refines the port neighbourhood (cost probe; weekly 2026-09-26, §10 chain F8)

**Status:** RUN `20261002T070008Z_ANS-4-step4b.log` — `[capture] rc=0`,
Status 0, 293 s, `memory.peak` 42 709 028 864 B (39.8 GiB), orphans none,
18 passed / 2 skipped (the skips are the record-fixture checks, which skip
by design off the 116 085-cell fixture, `:4865–4880`) (log footer; marked by
the 2026-10-02 03:00 review). **Price against the brackets (printed by the
review):** cells 259 543, INSIDE (250 k–270 k); unknowns ≈ 1 661 075, INSIDE
(1.6–1.8 M); four degree-2 drives 103.5 s, OUTSIDE, below the 250–450 s
bracket; 39.8 GiB, OUTSIDE, below the 80–160 GiB bracket (`:4843`,
`:4862–4863`). The case is cheaper than predicted on both axes. **Readout
(our numbers only):** the imported gates are green on both rungs. At degree 2
the reciprocity is 1.75e-13, the C4 spreads are 0.0285 / 0.0152 / 0.0188 %
against the 0.5 % band, and the column sums are 0.7173–0.7175 (`:4858–4861`).
Halving the conductor grading on the same cut moves the self-class \|S₁₁\| mean
by **+1.92 %**, the adjacent class by +0.47 % and the opposite class by
+0.61 %. The comparison is against entry 12's conductor-off `0.015:2` rung
(`20260927T070008Z_ANS-4-step3h.log:4793–4796`). On the same cut the global
ladder's port-1 \|S₁₁\| steps are +6.84 % (0.015 → 0.0075) and +3.50 %
(0.0075 → 0.005) (`…step3h.log:4793,6982,9257`). So the conductor move is
0.28× and 0.55× of those steps, in the same direction. The decision rule's
"comparable" and the commissioning of the `0.0075:2` / `0.005:2` rungs are
the 10-03 weekly's. The price came in below both brackets, which only makes
those rungs cheaper. *Was:* QUEUED 2026-09-28 for Friday 2026-10-02 02:00 —
`docs/testing/xl-queue.d/20261002-ANS-4-step4b.env` (18:00 daily review;
command copied verbatim). *Before that* READY (2026-09-28, `ANS-4` step 4 landed, §9 item 57) — the
knob name in the command below matches the test verbatim
(`FEM_EM_ANS4_STEP2_CONDUCTOR_FACTOR`), and the slot half's W2 meshed the
×0.5 congruent `0.015` rung at **259 543 cells**, inside the 250 k–270 k
bracket (`20260928T093258Z_ANS-4-step4-w2.log`); unset reproduced 2a″'s
digits to 8.2e-11 (`20260928T093139Z_ANS-4-step4-w1.log`). Was: PENDING
PREREQUISITE — `ANS-4` step 4 (chain F8's slot half, an
ordinary §9 item): `FEM_EM_ANS4_STEP2_CONDUCTOR_FACTOR` honoured on the
`RUNGSPEC` path (then read on the `RESOLUTION` path only), unset
digit-identical. `ANS-4` step 2 stays frozen; this is a new numbered step on
a knob that family never turned at degree 2.

**Licence class:** cost probe — no degree-2 point exists with the conductor
grading refined.

**Command** (entry 12's with the rungs shortened, the conductor factor added
and the raw-log name changed; final when the slot half lands):

```
XL_CHUNK="ANS-4-step4b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=10e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_CONDUCTOR_FACTOR=0.5 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step4b-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step4b-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step4b-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step4b-raw.log; cat /workspace/logs/ans4-step4b-raw.log; exit \$rc'"
```

**Price, predicted (brackets for the readout):** the halved conductor
grading on the congruent ×1 cut is measured at 259 543 cells (`PORT-21`
step 1b, degree 1); the nearest priced degree-2 rung is `ANS-4`'s
h = 0.0075 one — 281 728 cells, 1 819 094 unknowns, four drives 309.6 s
(3e) ⇒ **brackets: 250 k–270 k cells, 1.6–1.8 M unknowns, 250–450 s for
the degree-2 drives, 80–160 GiB**, printed INSIDE / OUTSIDE by the review.

**Readout (record):** cells, unknowns, times, `memory.peak`; the imported
gates on both rungs; the three classes at degree 2 beside entry 12's
`0.015:2` rung on the same cut with the conductor factor off — the move
under conductor refinement alone, at degree 2.

**Decision rule:** inside the brackets ⇒ the weekly commissions the same
knob on the `0.0075:2` and `0.005:2` rungs (the latter re-priced from this
point; it may be an `xxl` window); a self-class move under conductor
refinement comparable to the global ladder's steps ⇒ the conductor side is
the second unconverged axis and no `S₁₁` figure on a pinned-conductor
ladder is quoted as converged. A C4 spread outside the imported band on the
cut is a `GEO-32` finding for a §9 item, never a band change.

### 21. `xl` — `POST-7` step 2: `PORT-16`'s exact discrete power identity at degree 2 on F-human (weekly 2026-09-26, §10 chain H7)

**Status:** RUN `20261004T070008Z_POST-7-step2.log` — 13 passed, Status 0,
1097 s, 104.39 GiB; exact identity 1.207e-15 asserted at 1e-6, NEG 431×;
the identity-green branch fires (ledger row filled by the 2026-10-04 03:00
review). *Was:* QUEUED 2026-09-28 for Sunday 2026-10-04 02:00 —
`docs/testing/xl-queue.d/20261004-POST-7-step2.env` (18:00 daily review;
command copied verbatim; Saturday 10-03 is the `xxl` night). *Before that* READY (2026-09-28, `POST-7` step 1 landed, §9 item 58) — both
knob names in the command below match the test verbatim
(`FEM_EM_WF7_DEGREE`, `FEM_EM_WF7_IDENTITY_COLUMNS`). W2 is this command at
degree 1, `-n 8`: **279.7 s window (314 s pytest), 11.03 GiB summed
ru_maxrss**; (v) 0.0 on P17, P17 reproduction 2.4e-11 / 6.0e-8, NEG 481×
the band asserted (`20260928T094524Z_POST-7-step1-W2.log`); unset (W1) all 32
identity ≤ 2.6e-15, NEG 479× printed (`20260928T093759Z_POST-7-step1-W1.log`).
The route still runs the ≈ 125 s default-orientation mesh control after the
solve. Was: PENDING PREREQUISITE — `POST-7` step 1 (the degree knob on
`WF-7` step 1's gate module and the identity assert, degree-1 control
green in the heavy tier). **The knob names `FEM_EM_WF7_DEGREE` and
`FEM_EM_WF7_IDENTITY_COLUMNS` (how many driven columns are solved and
asserted; 1 = the first, one factorisation) on
`tests/validation/test_wf7_f_human_ring_matrix.py` are the contract**; if
step 1 lands other names this entry is re-issued, not edited.

**Why:** `xxl` 09-26 (`WF-7` step 0d) printed a terminal accounting
residual of 1.08–1.18e-2 on all 32 columns at degree 2 against
4.06–4.62e-3 at degree 1. `PORT-16` showed the terminal form misses by a
Cauchy–Schwarz deficit and that the statement which would convict an order
is a failure of the *exact discrete* identity — never yet asserted at
degree 2 on this fixture.

**Command** (final when step 1 lands; the module path is step 1's):

```
XL_CHUNK="POST-7-step2"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/post7-step2-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF7_DEGREE=2 FEM_EM_WF7_IDENTITY_COLUMNS=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, measured:** the same fixture and order, one drive — 3 256 418
unknowns, 322.95 s, 103.06 GiB `memory.peak` at `-n 16` (entry 9, 09-26);
the identity's assembly adds volume and sheet integrals ⇒ **predicted
450–900 s, 100–115 GiB**. Timeout 14 400 s; 512 GiB limit.

**Readout:** the exact identity's residual on the driven column **asserted** at the
imported `DISCRETE_IDENTITY_RTOL`; the terminal residual and the in-run κ
printed beside it; the negative control (one sheet's dissipation omitted)
asserted to miss by ≥ 10× the band.

**Decision rule:** identity green ⇒ the ≈ 1 % terminal residual is κ-class
at both orders and is registered as a record — degree 2 is not convicted;
identity red ⇒ known-issues entry, a degree-2 accounting diagnosis step, and
no human-scale degree-2 power figure is quoted until it closes. No band
moves.

### 22. `xl` — `PORT-21` step 3c: the κ-corrected sheet width at degree 2 on the h = 0.0075 congruent rung, 10 and 64 MHz (weekly 2026-09-30, §10 chain F12)

**Status:** SUPERSEDED 2026-10-03 (weekly) — its decision rule tests the
κ-corrected sheet width, a correction the operator's 2026-10-02 ruling no
longer applies; the exact-rectangle h-rung that answers T6 (i′) is entry 27.
Queue file removed in the weekly's commit; Monday 2026-10-05 is free.
*Was:* QUEUED 2026-10-01 (03:00 review) for Monday 2026-10-05 02:00 —
`docs/testing/xl-queue.d/20261005-PORT-21-step3c.env` (command line
byte-identical to the block below, checked with `cmp`). *Was:* READY
(2026-09-30) — every knob landed with `PORT-21` step 3
(F6, 2026-09-27; the three names below are exactly W3's). First free `xl`
night: Monday 2026-10-05 (the superseded entry 17's).

**Why:** tonight's private reading of F6's degree-2 window — at 10 MHz on
the ×1 cut, variant (a) takes the self-class miss against the order-matched
AED column from 15.0 % to 0.8 % — rests on one mesh. If the corrected value
is the converged-direction answer it must barely move under refinement
while the uncorrected one moves by its known step.

**Command** (entry 19's with the resolution, the chunk name and the raw-log
name changed — nothing else):

```
XL_CHUNK="PORT-21-step3c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/port21-step3c-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_PORT21_FREQUENCIES_MHZ=\"10 64\" FEM_EM_PORT21_DEGREE=2 FEM_EM_PORT21_RESOLUTION=0.0075 FEM_EM_PORT21_VARIANTS=\"control_prime a\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_port21_feed_sensitivity.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, measured in the family (scaling stated):** the `0.0075:2` rung of
the congruent ladder is 281 731 cells, 1 819 112 unknowns, four drives
280.0 / 303.9 s at `-n 16` (3h, 3i); F6's W3 measured three factorisations
per frequency (control′, the in-run κ solve, (a)) and a re-mesh per
frequency ⇒ **predicted 1 800–2 800 s**; memory by N^1.33 from the finest
rung's 280 GiB ⇒ **80–130 GiB** (one factor at a time, as W3 showed).
Timeout 14 400 s; 512 GiB limit.

**Readout (record; the module's imported gates asserted on every run):**
the in-run κ at each frequency; the diagonal of control′ and of (a) (the
`diagonal terms:` print, or the class means if `PORT-21` step 4 has landed
them in the JSON); (a)'s move against control′ per class. From them, by the
review's arithmetic against W3's ×1 values (10 MHz: control′
−0.3175+0.1375j, (a) −0.3799+0.1513j; 64 MHz from `PORT-21` step 4):
the successive change of the *corrected* and the *uncorrected* self class
between `0.015:2` and `0.0075:2`. Private: both variants against the
order-matched AED column.

**Decision rule (§10 T6, clause (i)):** the corrected sequence's successive
change on the self class ≤ ⅓ of the uncorrected one's on the same rungs
(the uncorrected steps are on record: 6.96 % at 10 MHz, 2.81 % at 64 MHz) ⇒
the κ-corrected width is the converged-direction convention and T6's clause
(ii) is read; otherwise tonight's 10 MHz reading was one mesh's coincidence,
it is written as one, and entries 17 / 18 are revived. **Predicted
(labelled, not a gate):** κ(0.0075:2) / κ(0.015:2) = 0.60 ± 0.10 at 10 MHz —
the ratio of the uncorrected misses on those rungs. No band moves; a gate
red is a `PORT-21` §9 item.

### 23. `xxl` — `PORT-21` step 3b: the κ-corrected sheet width at degree 2 on the finest congruent rung, 10 and 64 MHz (weekly 2026-09-30, §10 chain F7; entry 19 re-issued)

**Status:** RUN `20261003T070008Z_PORT-21-step3b.log` — Status 1, 9361 s,
291.36 GiB, both inside the price brackets; the one red is a missing
degree-1 h = 0.005 reference JSON, not a gate (known-issues OPEN 2026-10-03);
identity gates green at 10 / 64 MHz; κ printed; decision rule not fired (the
private half is the weekly's). Ledger row filled by the 2026-10-03 10:30
review. *Was:* QUEUED 2026-10-01 (03:00 review) for the `xxl` night of
Saturday 2026-10-03 02:00 — `docs/testing/xxl-queue.d/20261003-PORT-21-step3b.env`
(command line byte-identical to the block below, checked with `cmp`).
*Was:* READY (2026-09-30) — F6 landed the knob names exactly as below
(2026-09-27).

**Why this tier:** the job is `xl`-sized (≈ 280 GiB, ≈ 2.5 h). It is placed
on the `xxl` night because every `xl` night to 10-06 is taken, the 10-03
weekly reads chain F, and the `xxl` tier's own family (F-human degree 2)
has nothing ready — the window would otherwise go unspent a second week.

**Command** (entry 19's with the compose profile and service changed —
nothing else):

```
XL_CHUNK="PORT-21-step3b"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/port21-step3b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_PORT21_FREQUENCIES_MHZ=\"10 64\" FEM_EM_PORT21_DEGREE=2 FEM_EM_PORT21_RESOLUTION=0.005 FEM_EM_PORT21_VARIANTS=\"control_prime a\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_port21_feed_sensitivity.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Price, measured (F6's W3, 2026-09-27, scaled on the family's 1 476.9 s
sweep):** three factorisations per frequency and a re-mesh per frequency ⇒
**8 600–9 500 s, 275–295 GiB** (one factor at a time). Timeout 14 400 s
inside the 8 h window; 754 GiB limit.

**Readout:** as entry 22, on the finest rung; with entry 22 it makes the
corrected sequence three rungs long at 10 and 64 MHz.

**Decision rule, re-stated (entry 19's three branches assumed (a)'s
degree-2 move would be near its degree-1 value; at ×1 it is 5.4× it):**
read privately, (a)'s self class on this rung within 2 % of the
order-matched AED column at both frequencies **and** its change from the
×1 corrected value ≤ ⅓ of the uncorrected sequence's total change over the
same rungs ⇒ the sheet-width convention is the named, measured cause of the
self-class gap and `PORT-14`'s κ correction becomes the convention quoted
against AED (`PORT-21` step 5, no band moved); the corrected value drifting
with h as fast as the uncorrected one ⇒ the ×1 reading was a coincidence;
the corrected value stable but > 2 % from the column ⇒ κ is part of the gap
and F8's conductor axis is the rest.

### 24. `xl` — `ANS-7` step 4: Route B on the h = 0.005 hole mesh at degree 2, Cu column, the sparse grid (weekly 2026-09-30, §10 chain T12)

*(Superseded hold — the Status line below is current; this paragraph is history, kept verbatim. Marked by the 2026-10-04 03:00 review on the 10-03 weekly's flag 6.)* **HELD 2026-10-01 (10:30 review) — not issued.** A `log-pathologist`
ruling on item 132's logs (COUNTABLE) finds the degree-2 order move
κ-confounded: with κ held at its degree-1 value the raw degree-2 4-port S
carries ≈ 8 % unaccounted power (column sums ≈ 0.919, against 0.9986–0.9992
at degree 1), the size a degree-mismatched κ predicts, so the f₀ this
window would read against its 0.5 % decision rule carries an unmeasured κ
confound plausibly larger than 0.5 %. The command is issued only after T11b
(in-run degree-2 κ; serial on §9 item 140, F13) — then as a priced-family
variant of entry 14 with the in-run κ knob. Full ruling: PROJECT_PLAN §7
`ANS-7` row, 2026-10-01 10:30 annotation.

**Status: SUPERSEDED 2026-10-04 (10:30 review) by entry 28**, per the 10-03 weekly's rule (28 READY before the 10-07 18:00 review takes 24's night): queue file `20261008-ANS-7-step4.env` removed, `20261008-ANS-7-step4p.env` queued in its place. *Was:* QUEUED 2026-10-02 (10:30 review) for Thursday 2026-10-08 02:00 (`docs/testing/xl-queue.d/20261008-ANS-7-step4.env`, command copied verbatim). *Was:* **READY 2026-10-02 (06:00 CDT slot) — §9 item 173 landed green** (five-node Cu rows reproduce the seven-node rows to ≤ 1.5e-10 at p1 and ≤ 2.8e-10 at p2, against `BANDS` 1e-6; W2 / W3 solved = JSON ≤ 3.4e-11; `20261002T111807Z_ANS-7-item173-W2.log`, `…112024Z_ANS-7-item173-W3.log`). The command is unedited. *Was:* PENDING PREREQUISITE — §9 item 173; command issued below. The review read the module and found entry 24 not issuable as a knob variant. `test_assemble_from_node_json` assembles only on all seven default nodes and both columns, and `in_run` refuses `COLUMNS=cu` (`test_ans7_sparse_route_b.py:330–332,387,587–594`). At h = 0.005, degree 2, seven nodes × two columns plus the κ pass is ≈ 17 factorisations ≈ 22 000 s, over the 4 h cap. Item 173 adds an opt-in five-node, Cu-only, in-memory assembly (`FEM_EM_ANS7_SPARSE_ASSEMBLY=solved`), licensed at the default mesh against the seven-node rows. **The implementer landing 173 green marks this entry READY in the same commit and does not edit the command**; the next review queues it (first free `xl` night: Thu 10-08). The item-157 W0a flag was ratified by the 03:00 review (`docs/planning/chunks/ANS-7.md`).

**Command** (written 2026-10-02 03:00 review from the item-173 design; entry 14's durable-capture shape):

```
XL_CHUNK="ANS-7-step4"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans7-step4-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_SPARSE_MHZ=62.0,63.2,64.5,65.8,67.0 FEM_EM_ANS7_SPARSE_DEGREE=2 FEM_EM_ANS7_SPARSE_RESOLUTION=0.005 FEM_EM_ANS7_SPARSE_COLUMNS=cu FEM_EM_ANS7_SPARSE_FIELD_SOLVE=0 FEM_EM_ANS7_SPARSE_KAPPA=in_run FEM_EM_ANS7_SPARSE_ASSEMBLY=solved timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans7_sparse_route_b.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the run writes only under the gitignored `/workspace/logs/` (`ans7_sparse/kappa_p2_h0.005.json` and five `f_*_p2_h0.005_kin_cu.json`).

**Price, measured in the family** (entry 14, `20260929T070008Z_ANS-6-step2c.log`, `-n 16`):
- **κ pass:** entry 14's own route, ≈ 4 105 s. That is a mesh (104.6 s), the PEC factorisation (1 461.2 s), the Cu factorisation (1 302.3 s) and the P1 field solve (1 236.9 s). `FIELD_SOLVE=0` cannot shorten it, because κ exists only in the field-solve branch.
- **Grid:** a second mesh (≈ 105 s) plus five Cu sweeps at 1 237–1 461 s each.
- **Total ⇒ predicted 10 400–11 600 s, 260–300 GiB**, one factor at a time. Item 173's `del klad; gc.collect()` keeps the κ factor from overlapping the grid's; entry 14 peaked at 266.79 GiB. Timeout 14 400 s, limit 512 GiB.

**Readout (record):**
- f₀, Q, R_in and the k = ±1 share of the Cu column, from the interpolated, terminated Z(f) on five nodes.
- Beside them: the ×1 degree-1 rows, the ×1 degree-2 in-run rows (62.704216 MHz, 18 570 Ω, 984.2) and the h-move.
- κ at p2, h = 0.005, pooled and per sheet; (R2) against entry 14's line; (N).
- Δ_LOO and B₅; Im Y_in at every node; the 3-node control.
- Private: the rows against the AED columns.

**Decision rule:** as written below. Additions for the node set:
- A reading within B₅ of the 0.5 % edge is undecided.
- Crossings ≠ 1 (Im Y_in(62.0) > 0) ⇒ f₀ lies below the window. That bound is the finding; the review re-poses the grid and the `_kin_cu` JSONs are reused.
- (R2) or (N) red ⇒ κ is unlicensed on this mesh. The rows are recorded as unlicensed and the reading becomes a §9 item; no band moves.

**Predicted (labelled, not a gate):** f₀,Cu between 62.70 and 63.98 MHz.

*Was:* **PREREQUISITE LANDED 2026-10-02 — knob `FEM_EM_ANS7_SPARSE_KAPPA=in_run` exists; command to be issued by the next review** (§9 item 157, T11b; κ_deg2 = 5.2094e-2 on the default mesh, reading "κ confound"). κ is mesh-dependent, so this entry runs `in_run`, never `stored`, and its κ pass includes the PEC sweep. **Flag for the review:** W0a's knob-unset print missed digit identity (see the item 157 outcome); confirm it before issuing.

*Was:* **PREREQUISITE LANDED 2026-10-01 — knob names `FEM_EM_ANS7_SPARSE_*` exist
(`_MHZ`, `_DEGREE`, `_RESOLUTION`, `_FIELD_SOLVE`, `_COLUMNS`, in
`tests/validation/test_ans7_sparse_route_b.py`); command to be issued by the
next review.** (§9 item 132; degree-2 Δ_LOO ≤ 8.2e-11 inside the 1e-6
bands on the 0.015 m mesh; the order move R_in / Q ≈ −98 % is printed and
uninterpreted — the review should read it before pricing.)

**Status (as written):** PENDING PREREQUISITE — `ANS-7` step 3 (chain T11, an ordinary
§9 item): the sparse-grid Route B with its degree, resolution and
frequency-grid knobs, the degree-1 identity against the 51-point sweep
green and the ×1 degree-2 order move printed. **The command is issued when
those knob names exist** (by the weekly, or by the daily as a priced-family
variant of entry 14 — same fixture and mesh family, the frequency list
varied); it must carry `XL_OUTPUTS` if it writes into the tree.

**Why:** `ANS-7`'s absolute rows were ruled INCONCLUSIVE on one degree-1
rung (`f₀` 1.2–1.9 % above AED). This is the order- and h-matched rung.

**Price, measured:** entry 14 — 523 816 cells, 3 466 952 unknowns, the Cu
factorisation 1 302 s, 266.79 GiB at `-n 16`. A four-port `Z` at five
frequencies is five factorisations ⇒ **predicted 6 500–8 000 s,
260–290 GiB**. Timeout 14 400 s.

**Readout:** `f₀`, Q, `R_in` and the k = ±1 share from the interpolated,
terminated `Z(f)`, beside the ×1 degree-1 and degree-2 values. Private: the
three rows against the AED columns.

**Decision rule:** `f₀` on this rung within 0.5 % of the order-matched
First Order column (AED's own order move) ⇒ the 09-30 finding was
discretisation and `ANS-7` is ruled AGREE at the stated level; a gap that
does not close with order and h ⇒ a model difference, read with chain F's
κ result before any diagnosis chunk is opened.

### 25. `xl` — `ANS-4` step 4c: entry 20's conductor-refined degree-2 rung at 64 MHz (daily licence: priced-family variant, 2026-10-02 03:00 review; §10 chain F8)

**Status:** RUN `20261006T070008Z_ANS-4-step4c.log`: 18 passed, 2 skipped, Status 0, 318 s, 40.75 GiB, 259 543 cells; self-class \|S₁₁\| conductor move +0.48 % at 64 MHz against 3i's global step +2.35 % (ledger row, 2026-10-06 03:00 review). *Was:* QUEUED 2026-10-02 (03:00 review) for Tuesday 2026-10-06 02:00 —
`docs/testing/xl-queue.d/20261006-ANS-4-step4c.env` (command copied
verbatim). *Before that* READY on writing: entry 20's knobs all exist on
`main`, so there is no prerequisite.

**Licence class:** priced-family variant of entry 20. It uses the same
module, mesh, `RUNGSPEC`, conductor factor and rank count, and varies **only
the frequency**: `FEM_EM_ANS4_FREQUENCY_HZ=64e6`, in place of 10e6. Its
partner without the conductor factor is entry 13 (3i), the same cut at
64 MHz. **Why:** F8 asks whether the conductor grading is a second
unconverged axis behind the `S₁₁` gap. That gap is open at 10 **and** 64 MHz
(§2), and entry 20 answers only at 10 MHz.

**Command** (entry 20's with the frequency, the chunk name and the raw-log
name changed — nothing else):

```
XL_CHUNK="ANS-4-step4c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=64e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_CONDUCTOR_FACTOR=0.5 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step4c-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step4c-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step4c-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step4c-raw.log; cat /workspace/logs/ans4-step4c-raw.log; exit \$rc'"
```

**Price, measured in the family (scaling stated):** entry 20 measured 293 s
end to end, four degree-2 drives in 103.5 s and 39.8 GiB at `-n 16`
(`20261002T070008Z_ANS-4-step4b.log:4862–4863`). Only the frequency changes,
so the mesh and the unknowns are identical. On the congruent ladder a direct
factorisation's price barely moves with frequency: the `0.0075:2` rung's four
drives took 280.0 s at 10 MHz (3h) and 303.9 s at 64 MHz (3i), as entry 22
records. **Brackets: 259 543 cells
(asserted equal by the rung's own print), 200–600 s, 30–60 GiB.**

**Readout (record):** the imported gates on both rungs. Then the three
classes at degree 2 beside entry 13's `0.015:2` rung at 64 MHz
(`20260928T070008Z_ANS-4-step3i.log`), which is the move under conductor
refinement alone. The review prints that move beside 3i's global-ladder
steps on the same cut, as entry 20's readout did at 10 MHz.

**Decision rule:** entry 20's, read at 64 MHz. If the self-class move under
conductor halving is comparable to the global ladder's steps at 64 MHz, the
conductor side is a second unconverged axis at both open-gap frequencies,
and no 64 MHz `S₁₁` figure on a pinned-conductor ladder is quoted as
converged. A C4 spread outside the imported band is a `GEO-32` finding for
a §9 item, never a band change. The rung extension (`0.0075:2` / `0.005:2`)
stays the weekly's to commission.

### 26. `xl` — `ANS-4` step 4d: entry 20's conductor-refined degree-2 rung at 128 MHz (daily licence: priced-family variant, 2026-10-02 03:00 review; §10 chain F8)

**Status:** RUN `20261007T070008Z_ANS-4-step4d.log`: 18 passed, 2 skipped, Status 0, 306 s, 40.25 GiB, 259 543 cells; self-class \|S₁₁\| conductor move +0.24 % at 128 MHz against 3j's global step +1.02 % (ledger row, 2026-10-07 03:00 review). *Was:* QUEUED 2026-10-02 (03:00 review) for Wednesday 2026-10-07 02:00 —
`docs/testing/xl-queue.d/20261007-ANS-4-step4d.env` (command copied
verbatim). *Before that* READY on writing (no prerequisite).

**Licence class:** priced-family variant of entry 20, with the frequency
only varied: `FEM_EM_ANS4_FREQUENCY_HZ=128e6`. Its partner without the
conductor factor is entry 15 (3j). **Why:** this completes F8's conductor
reading at the third frequency `ANS-4`'s ladder ran. At 128 MHz the
order-matched couplings and `S₁₁` already agree (§2), so a large conductor
move here would say the agreement there is not converged either.

**Command** (entry 20's with the frequency, the chunk name and the raw-log
name changed — nothing else):

```
XL_CHUNK="ANS-4-step4d"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_FREQUENCY_HZ=128e6 FEM_EM_ANS4_STEP2_C4_CONGRUENT=1 FEM_EM_ANS4_STEP2_CONDUCTOR_FACTOR=0.5 FEM_EM_ANS4_STEP2_RUNGSPEC=\"0.015:1 0.015:2\" timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_resolution_ladder.py -v -s --tb=short > /workspace/logs/ans4-step4d-raw.log 2>&1; rc=\$?; echo \"[XL] memory.peak bytes:\" >> /workspace/logs/ans4-step4d-raw.log; cat /sys/fs/cgroup/memory.peak >> /workspace/logs/ans4-step4d-raw.log; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans4-step4d-raw.log; cat /workspace/logs/ans4-step4d-raw.log; exit \$rc'"
```

**Price, measured in the family:** as entry 25. **Brackets: 259 543 cells,
200–600 s, 30–60 GiB.**

**Readout (record):** as entry 25, beside entry 15's `0.015:2` rung at
128 MHz (`20260930T070008Z_ANS-4-step3j.log`) and 3j's global steps.

**Decision rule:** as entry 25, read at 128 MHz.

### 27. `xl` — `ANS-4` step 6b: the h = 0.0075 rung on exact SPEC rectangles at degree 2, κ off, 10 / 64 / 128 MHz (weekly 2026-10-03, §10 chain F18b)

**Status:** RUN `20261005T070009Z_ANS-4-step6b.log`: 16 passed and 5 skipped, Status 0, 1 635 s, 90.76 GiB, 282 198 cells. (G) is green κ off at all three frequencies. Self class move against ×1 is 0.25 / 0.23 / 0.22 % at 10 / 64 / 128 MHz (`:3726–3728`). The 2026-10-05 03:00 review filled the ledger row; the private column clause is the weekly's. *Was:* QUEUED 2026-10-04 (10:30 review) for Monday 2026-10-05 02:00 —
`docs/testing/xl-queue.d/20261005-ANS-4-step6b.env` (command issued below by
this review and copied verbatim). *Before that* READY (2026-10-04, §9 item 212): knobs present, unset
digit-identical (W0 log `20261004T093239Z_ANS-4-item212-W0-p1-unset.log`,
21 passed, (R) ≤ 9.6e-11), resolution knob reaches the mesh (W0b log
`20261004T093512Z_ANS-4-item212-W0b-p1-r0p0075.log`, 282 198 cells at
h = 0.0075 against 117 206 unset). The contract was the knobs
`FEM_EM_ANS4_EXACT_RESOLUTION` and `FEM_EM_ANS4_EXACT_CONDUCTOR_FACTOR` on
the exact-strip module (`tests/validation/test_ans4_exact_strip_records.py`,
beside the existing `FEM_EM_ANS4_EXACT_DEGREE`), unset digit-identical to
item 192's records.
**Those names are the contract; the command is issued when they exist** (by
the daily, as entry 20's durable-capture shape with the module, the knobs
and the raw-log name changed), with `FEM_EM_ANS4_EXACT_RESOLUTION=0.0075`,
`FEM_EM_ANS4_EXACT_DEGREE=2`, all three frequencies, `-n 16`.

**Why:** tonight's private reading puts ours at degree 2 on exact rectangles
within 0.73 % of AED First Order on every class at 10 / 64 / 128 MHz — on
one mesh. T6 clause (i′) needs the same figure one rung finer.

**Price, measured in the family (scaling stated):** the `0.0075:2` rung of
the congruent ladder is 281 731 cells, 1 819 112 unknowns, four drives in
280.0 / 303.9 s at `-n 16` (3h, 3i); three frequencies and a mesh ⇒
**predicted 1 200–2 200 s, 80–130 GiB** (N^1.33 from the finest rung's
280 GiB). The exact-strip mesh differs from the congruent cut only at the
port boxes. Timeout 14 400 s; 512 GiB limit.

**Readout (record; `PORT-11`'s imported gates asserted on the rung):** the
class means at degree 2 per frequency beside item 192's ×1 values
(`20261003T094036Z_ANS-4-item192-W3-p2.log:2199–2216`); the successive
change per class; in-run κ printed as a diagnostic. Private: the rung
against the First Order column.

**Decision rule (§10 T6, clause (i′)):** the self class's successive change
≤ 2.32 % at 10 MHz and ≤ 0.94 % at 64 MHz (⅓ of the ragged uncorrected steps
on the same rungs, the 09-30 thresholds) **and** the rung still inside 2 %
of the column ⇒ degree 2 on exact rectangles is the licensed convention and
`PORT-21` closes; a larger change ⇒ the ×1 agreement was one mesh's and the
h = 0.005 rung is priced from this one. No band moves.

**Command** (issued 2026-10-04 10:30 review: entry 24's durable-capture shape
— orphan counts before and after — with the module, the knobs and the raw-log
name changed; `PYTHONPATH` as item 212's windows, because the module imports
`tests.*`; `FEM_EM_ANS4_EXACT_MHZ` unset ⇒ all three frequencies on one mesh,
κ in-run at 128 MHz printed as the diagnostic, `test_ans4_exact_strip_records.py:434–446`;
conductor factor unset ⇒ the ×1 conductor grading):

```
XL_CHUNK="ANS-4-step6b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans4-step6b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_EXACT_DEGREE=2 FEM_EM_ANS4_EXACT_RESOLUTION=0.0075 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_exact_strip_records.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the knobbed path writes only the gitignored
`logs/ans4_exact/ans4_exact_p2_r0.0075_cx_f10-64-128.json`-style file
(`test_ans4_exact_strip_records.py:560`), never item 192's.

### 28. `xl` — `ANS-7` step 4′: Route B on the h = 0.005 hole mesh with exact rectangles, κ off, degree 2, Cu, five nodes (weekly 2026-10-03, §10 chain T12′; replaces entry 24)

**Status:** RUN `20261008T070008Z_ANS-7-step4p.log` — 17 passed, 18 skipped,
Status 0, 7 675 s (2.3 % above the bracket), 272.35 GiB, 523 880 cells; one
up-crossing; Cu f₀ 62.751784 MHz, h-move against ×1 f₀ +0.0809 %, R_in
−1.185 %, Q −1.199 % (ledger row; marked 2026-10-08 10:30 review). The
decision rule's private comparison is the 10-10 weekly's. *Was:* QUEUED
2026-10-04 (10:30 review) for Thursday 2026-10-08 02:00 —
`docs/testing/xl-queue.d/20261008-ANS-7-step4p.env` (command issued below by
this review and copied verbatim); entry 24 is superseded and its queue file
removed. *Was:* READY (2026-10-04, §9 item 214, `ANS-7` step 6, chain T19
green) — command not yet issued; the next review issues it (entry 24's
durable-capture shape). Knobs on `tests/validation/test_ans7_exact_strip.py`:
`FEM_EM_ANS7_EXACT_ROUTE=B`, `FEM_EM_ANS7_EXACT_DEGREE=2`,
`FEM_EM_ANS7_EXACT_RESOLUTION=0.005`, `FEM_EM_ANS7_EXACT_COLUMNS=cu`,
`FEM_EM_ANS7_EXACT_MHZ=62.0,63.2,64.5,65.8,67.0` (unset resolution reproduces
item 196's rows to ≤ 2.2e-12, (R2); the Cu-only five-node set is licensed at
×1, (S5): rows vs seven-node ≤ 9.2e-10, Δ_LOO ≤ 3.5e-8). Node JSONs land as
`logs/ans7_exact/node_<mhz>_p2_h0.005_cu.json`; the assembly runs in the
window that completes the five. **Measured per-node price on the T19 rung**
(h = 0.01, 140 270 cells, 968 016 degree-2 unknowns, `-n 12`): Cu sweep
143.5–159.8 s per node, mesh 53.6–55.4 s (built twice per window: tie +
ladder), summed `ru_maxrss` 32.8–36.5 GiB
(`20261004T111415Z_ANS-7-item214-W1-probe.log`, `…112040Z_…-W2-rung.log`,
`…112945Z_…-W3-rung.log`). **If this entry is READY by the 2026-10-07 18:00
review it takes entry 24's night (2026-10-08) and entry 24 is superseded;
otherwise entry 24 runs as queued** — its ×1 twin reads the same `f₀` to
0.005 %, so its h-move is admissible.

**Price, measured in the family:** entry 14's Cu factorisation on this mesh,
1 302 s at `-n 16`, 266.79 GiB; five nodes and a mesh, no κ pass (κ is off)
⇒ **predicted 6 300–7 500 s, 260–300 GiB** — ≈ 4 100 s cheaper than entry
24, which pays for an in-run κ. Timeout 14 400 s; 512 GiB limit.

**Readout (record):** `f₀`, Q, `R_in` and the k = ±1 share of the Cu column
beside item 196's ×1 rows (62.701082 MHz / 18 441.75 Ω / 977.26) and T19's
intermediate rung; the h-move; Δ_LOO; Im `Y_in` at every node. Private: the
rows against the First Order column.

**Decision rule (entry 24's, unchanged):** `f₀` within 0.5 % of the
order-matched First Order column on this rung ⇒ the 09-30 finding was
discretisation and `ANS-7` is ruled AGREE at the stated level, with the
h-move quoted beside it; a gap that opens with h ⇒ the ×1 agreement was a
cancellation and is written as one. Crossings ≠ 1 ⇒ the bound is the
finding and the grid is re-posed.

**Command** (issued 2026-10-04 10:30 review: entry 24's durable-capture shape
with the module and knobs of item 214's windows, the resolution set to 0.005,
all five Cu nodes in one window so the assembly runs, and the raw-log name
changed; `FEM_EM_ANS7_EXACT_TUNE` unset, so item 217's T21 tests skip):

```
XL_CHUNK="ANS-7-step4p"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans7-step4p-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_ROUTE=B FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_RESOLUTION=0.005 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_MHZ=62.0,63.2,64.5,65.8,67.0 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: node JSONs land under the gitignored
`logs/ans7_exact/node_<mhz>_p2_h0.005_cu.json` (`test_ans7_exact_strip.py:195–203`).

### 29. `xl` — `WF-7` step 0e: the F-human σ = 800 32×32 at degree 2, κ off, saved (weekly 2026-10-03, §10 chain H15b)

**Status:** RUN `20261009T070009Z_WF-7-step0e-H15b.log`: 16 passed, 1 skipped (split control, by design), Status 0, 782 s, 103.95 GiB, 507 266 cells; reciprocity, σ_max and C16 classes green at degree 2; matrix saved (ledger row, 2026-10-09 18:00 review). The decision rule is the 10-10 weekly's (private comparison). *Was:* QUEUED 2026-10-05 (03:00 review) for Friday 2026-10-09 02:00 —
`docs/testing/xl-queue.d/20261009-WF-7-step0e-H15b.env` (command issued below
by this review and copied verbatim). *Was:* READY (2026-10-04 19:30 slot, §9 item 231; prerequisite §10
H15a landed) — the knob `FEM_EM_WF7_SAVE_NPZ` on
`tests/validation/test_wf7_f_human_ring_matrix.py` writes the held 32×32
under the gitignored `/workspace/logs/`, and the degree-1 matrix reproduces
through the reloaded file (σ_max and all 18 class means within
`REUSE_REPRODUCTION_RTOL`, `20261005T003402Z_WF-7-H15a-r3-W2.log`). Knob
convention: `FEM_EM_WF7_ORIENTATION_CONTROL` unset is the landed inline path
(control built after the matrix), `split` moves it out of the matrix window,
`only` runs it alone. **The command is still the review's to issue** (entry
21's shape with `FEM_EM_WF7_DEGREE=2`, all 32 columns, the save knob set).

**Why:** `ANS-8` rows 1 / 4 are INCONCLUSIVE at degree 1 — our `f₀` sits
1.3–1.9 % above the σ = 800 columns, the signature that on F-small turned
out to be element order. `xxl` 09-26 solved this matrix and did not keep it.

**Price, measured:** entry 9 (`WF-7` 0d, `xxl` 2026-09-26) — 646 s, 103 GiB
at `-n 16` for all 32 columns at degree 2 ⇒ **predicted 600–900 s,
100–115 GiB**: an `xl`-sized job. Timeout 14 400 s; 512 GiB limit.

**Readout (record; the module's imported gates asserted):** reciprocity,
σ_max, the C16 class spreads at degree 2; `C_mode`, C(k = 2) and the
feed-ring share / leak through `ring_mode_rows_v2`, beside the degree-1 rows
(`20261002T050705Z_ANS-8-v2-W1-R.log`); the order move. Private: rows 1 / 4
against the First Order column.

**Decision rule:** the degree-2 k = 1 resonance within 0.5 % of the First
Order σ = 800 column ⇒ `ANS-8` rows 1 / 4 are ruled AGREE at degree 2 and
human-scale tuned figures are quoted at degree 2 only; a gap that survives
the order ⇒ a human-scale finding, known-issues, read with H14's sheet
census before any diagnosis. **Predicted (labelled, not a gate):** `f₀`
falls 1.5–2.5 % from degree 1.

**Command** (issued 2026-10-05 03:00 review). It uses entry 21's durable-capture shape with
`FEM_EM_WF7_DEGREE=2`, `FEM_EM_WF7_IDENTITY_COLUMNS` unset (all 32 columns;
the knob is exclusive with the save knob, `test_wf7_f_human_ring_matrix.py:221–223`),
`FEM_EM_WF7_SAVE_NPZ` under the gitignored `/workspace/logs/wf7/`, and
`FEM_EM_WF7_ORIENTATION_CONTROL=split`. That is the configuration item 231's W2
ran green at degree 1, `-n 8` (`20261005T003402Z_WF-7-H15a-r3-W2.log:12`).
The degree-1 records (P17 reproduction, reloaded σ_max and class means) are
printed only at degree 2 (`:634–635`, `:719`, `:841`, `:867`). The
`WINDOW_BUDGET_S` 540 s stop rule (`:583`) restricts the identity to drives
1 / 9 / 17 / 25 if the projected window exceeds it. That is a printed fallback,
not a red. κ is not applied (`kappa_applied=False`, `:670`).

```
XL_CHUNK="WF-7-step0e-H15b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs/wf7 && R=/workspace/logs/wf7-h15b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF7_DEGREE=2 FEM_EM_WF7_ORIENTATION_CONTROL=split FEM_EM_WF7_SAVE_NPZ=/workspace/logs/wf7/wf7-s32-degree2-n16.npz timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the saved 32×32 lands under the gitignored
`/workspace/logs/wf7/wf7-s32-degree2-n16.npz`; H12c / the `ANS-8` re-read
load it from there.

### 30. `xl` — `ANS-4` step 6c: the h = 0.005 rung on exact SPEC rectangles at degree 2, κ off, 10 / 64 / 128 MHz (weekly 2026-10-05, §10 chain F25)

**Status:** READ 2026-10-12 (weekly, §10 F27) — the decision rule's first branch fired: from the 10-digit S lines the self-class 0.0075 → 0.005 change is **0.0144 / 0.0569 / 0.0662 %** at 10 / 64 / 128 MHz (adjacent 0.0088 / 0.0235 / 0.0472 %, opposite 0.0161 / 0.0611 / 0.1110 %; ratios to the ×1 → 0.0075 move 0.057 / 0.243 / 0.301 against the p = 2 prediction 0.185), all self-class changes ≤ 0.10 % ⇒ §2 quotes the degree-2 exact figure with the h-statement "h-converged to ≤ 0.1 % on three rungs"; privately the rung stays within 1 % of AED First Order on every class (the 10 MHz self class is the worst and is flat with h), so the public level is unchanged and no known-issues entry is opened (`docs/private/f27-f26-h15c-ans9-adjudication-2026-10-12.md`). *Was:* RUN `20261011T070008Z_ANS-4-step6c.log`: 16 passed, 5 skipped, Status 0, 7 347 s (22 % above the bracket), 288.09 GiB, 593 478 cells; (G) green at all three frequencies; self-class move against ×1 0.2519 / 0.2868 / 0.2866 % at 10 / 64 / 128 MHz — F27 is the 10-14 weekly's (ledger row, 2026-10-11 10:30 review). *Was:* QUEUED 2026-10-06 (03:00 review) for Sunday 2026-10-11 02:00 — `docs/testing/xl-queue.d/20261011-ANS-4-step6c.env` (command copied verbatim). *Was:* READY 2026-10-05 (weekly) — the knobs exist (`FEM_EM_ANS4_EXACT_RESOLUTION`, `FEM_EM_ANS4_EXACT_DEGREE`, `test_ans4_exact_strip_records.py:114–117`; entry 27 ran them at 0.0075 on 10-05) and the command below is entry 27's with the resolution and the raw-log name changed. The clerk queues it for the first free night: **Sunday 2026-10-11 02:00**.

**Why:** T6 is decided on two h-rungs (×1, h = 0.0075) and one conductor rung; what is **not** licensed is an h-converged degree-2 figure. On the self class the ×1 → 0.0075 move is 0.25 / 0.23 / 0.22 % at 10 / 64 / 128 MHz, and at 10 MHz it is *away* from the private column (§10 rulings 2026-10-05). A third rung makes the convergence statement §2 can carry, or shows there is none.

**Price, measured in the family (scaling stated):** the ragged `0.005:2` rung of the congruent ladder is 592 550 cells, 3 827 268 unknowns, four drives in 1 378.6 / 1 563.1 s at `-n 16`, ≈ 275–282 GiB `memory.peak` (3h, 3i, 3j); the exact-strip mesh differs only at the port boxes (282 198 against 281 731 cells at 0.0075). Three frequencies on one mesh plus the mesh ⇒ **predicted 4 500–6 000 s, 260–300 GiB**. Timeout 14 400 s; 512 GiB limit.

**Readout (record; `PORT-11`'s imported gates asserted on the rung):** the class means at degree 2 per frequency beside item 192's ×1 values and entry 27's h = 0.0075 values; the successive change per class; in-run κ as a diagnostic. Private: the rung against the First Order column.

**Decision rule (pre-registered; §10 F27, read by the 10-14 weekly).** Step ratios 0.015 : 0.0075 : 0.005 on the global `resolution` knob (the conductor grading is pinned, as in entry 27). For a second-order method the 0.0075 → 0.005 self-class move is (0.0075² − 0.005²) / (0.015² − 0.0075²) = 0.185 × the ×1 → 0.0075 move, i.e. ≈ 0.05 %. **Observed ≤ 0.10 % (twice the p = 2 prediction) on the self class at every frequency ⇒ §2 quotes the degree-2 exact figure with the h-statement "h-converged to ≤ 0.1 % on three rungs"; 0.10–0.25 % ⇒ the h-statement is the measured bracket and the public wording stays "within 1 % of AED First Order on three meshes"; > 0.25 % or non-monotone ⇒ a known-issues finding and a diagnosis step.** No band moves either way.

**Command** (entry 27's durable-capture shape with `FEM_EM_ANS4_EXACT_RESOLUTION=0.005` and the raw-log name changed; `FEM_EM_ANS4_EXACT_MHZ` unset ⇒ all three frequencies on one mesh; conductor factor unset ⇒ the ×1 conductor grading):

```
XL_CHUNK="ANS-4-step6c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans4-step6c-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_ANS4_EXACT_DEGREE=2 FEM_EM_ANS4_EXACT_RESOLUTION=0.005 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans4_exact_strip_records.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the knobbed path writes only the gitignored
`logs/ans4_exact/ans4_exact_p2_r0.005_cx_f10-64-128.json`-style file, never
item 192's records.

### 31. `xl` — `ANS-6` step 5b: the exact hole route at h = 0.0075, degree 2, Cu and PEC, 64 MHz (weekly 2026-10-05, §10 chain F26)

**Status:** READ 2026-10-12 (weekly, §10 F26 / F27) — both clauses of the decision rule met: the Cu − PEC class gap on the rung differs from F19's ×1 value by ≈ 3 % of itself (0.0615 vs 0.0633 %, 0.0279 vs 0.0287 %, 0.0479 vs 0.0495 % — inside ⅓), and privately the Cu and PEC rungs sit within 0.21 % of the AED Finite-Conductivity / PEC First Order columns on every class at 64 MHz (×1: 0.33 %; and F19's ×1 columns at 10 / 128 MHz, read tonight for the first time, within 0.10 / 0.31 %) ⇒ `ANS-6` row (a) is written as an h-statement at 64 MHz and the row goes ✅ (its glyph was held for this re-read); entry 34 (F31) stays queued as the copper-loss partition's third rung for T30 (`docs/private/f27-f26-h15c-ans9-adjudication-2026-10-12.md`). *Was:* RUN `20261012T070009Z_ANS-6-step5b.log` — 25 passed, 14 skipped, Status 0, 1 532 s, 76.17 GiB, 226 904 cells; copper identity green; readout in the ledger row (2026-10-12 03:00 review); decision rule and F27 are the weekly's. *Was:* QUEUED 2026-10-06 (18:00 review) for Monday 2026-10-12 02:00 — `docs/testing/xl-queue.d/20261012-ANS-6-step5b.env` (command below, issued by this review and copied verbatim). *Was:* READY (2026-10-06, §9 item 266): knob present, unset digit-identical (W0 / W1 / W2 logs), resolution knob reaches the mesh (WK log, cells). *Was:* PENDING PREREQUISITE — §10 F24 (`ANS-6` step 5b′): a resolution knob `FEM_EM_EXACT_ANS6_RESOLUTION` on `tests/validation/test_ans6_exact_strip.py`, beside the existing `FEM_EM_EXACT_ANS6_PART` / `_DEGREE` / `_MODE`, unset digit-identical to F19's records (items 226 / 232) and shown to reach the mesh at 0.0075. **The implementer landing F24 green marks this entry READY in the same commit and does not write the command**; the next review issues it in entry 27's durable-capture shape with the module, the knobs (`FEM_EM_EXACT_ANS6_RESOLUTION=0.0075`, `FEM_EM_EXACT_ANS6_DEGREE=2`, the module's frequency knob at 64 only) and the raw-log name changed, `-n 16`, and queues it for **Monday 2026-10-12 02:00** or the first free night after.

**Why:** F19 recorded the exact hole route at degree 2 on one mesh (×1) and the copper column was read privately on it (row (a)); the conductor's part of the tuned 1/Q sits 4.4 % high at both orders on `ANS-7`. One h-rung on the copper hole mesh says whether the copper level is an h-statement at 64 MHz or still one mesh's.

**Price, measured in the family (scaling stated):** entry 16 (`ANS-6` 2d, the ragged hole mesh at h = 0.0075, degree 2, 64 MHz, PEC / Cu / P1) ran 818 s; entry 14 (2c, h = 0.005) 4 110 s at 266.79 GiB with 523 816 cells. The 0.0075 mesh is ≈ 0.3× the cells ⇒ N^1.33 gives ≈ 55 GiB. **Predicted 800–1 500 s, 40–90 GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; `TH-14`'s imported gates and the copper surface-loss identity asserted on the rung):** Cu and PEC class means at degree 2 beside F19's ×1 values; the Cu − PEC move per class on each rung; the copper surface-loss partition on the rung. Private: the rung against the Finite-Conductivity and PEC columns.

**Decision rule (pre-registered; §10 F27, read by the 10-14 weekly):** the Cu − PEC class moves on the h = 0.0075 rung within ⅓ of F19's ×1 values of the same moves **and** the rung inside the 2 % column clause ⇒ `ANS-6`'s row (a) is written as an h-statement at 64 MHz; larger ⇒ the h = 0.005 rung is priced from this one (entry 14's 4 110 s says it fits). No band moves.

**Command** (issued 2026-10-06 18:00 review: entry 30's durable-capture shape with item 266 W1's module and knobs — `FEM_EM_EXACT_ANS6=solve`, `FEM_EM_EXACT_ANS6_DEGREE=2`, `TH14_STEP2_FREQS_MHZ=64` (`20261006T183450Z_ANS-6-item266-W1-p2-64-identity.log:12`) — plus `FEM_EM_EXACT_ANS6_RESOLUTION=0.0075`; `FEM_EM_EXACT_ANS6_PART` unset ⇒ the item-193 / item-216 call, sweeps and identity on one mesh, so the Cu and PEC class means and the copper identity both land; the knobbed path asserts (F′) and skips (R) / (N3) by design):

```
XL_CHUNK="ANS-6-step5b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans6-step5b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 FEM_EM_EXACT_ANS6_RESOLUTION=0.0075 TH14_STEP2_FREQS_MHZ=64 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the knobbed path writes only the gitignored `_r0.0075` JSON, never F19's records.

### 32. `xl` — `WF-12` step 4: one degree-2 copper column on the F-human hole mesh, κ off, 64 MHz, `-n 16` (weekly 2026-10-07, §10 chain H16)

**Status:** RUN `20261013T070008Z_WF-12-step4-H16.log` — 13 passed, Status 0, 514 s, 89.31 GiB `memory.peak`; build 122.0 s, factorising drive 287.1 s (solve 273.1 + post 14.0), held re-solve 3.76 s; (P) 1.468e-13, NEG 119.0×, (Col) 0.9849 green; S17,17 d1 → d2 move 3.550 %, coil share 0.3775 (d1 0.3704); the degree-2 column is not saved to disk (ledger row, 2026-10-13 03:00 review). **Decision rule's first branch on its numbers:** identities green and the factorising drive 287 s ≤ 1 500 s, memory ≤ 200 GiB ⇒ the degree-2 copper 32×32 is an `xl` case priced ≈ 122 + 287 + 31 × 5.64 ≈ 584 s (`:6762`) — entry 36 is re-tiered to `xl` by the 10-14 weekly per its own rule. *Was:* QUEUED 2026-10-08 (10:30 review) for Tuesday 2026-10-13 02:00 — `docs/testing/xl-queue.d/20261013-WF-12-step4-H16.env` (command issued below by this review and copied verbatim; trailing-7-day count on 10-13 is six with it). *Was:* READY 2026-10-08 (§9 item 312). H16a landed: `tests/validation/test_wf12_f_human_copper_degree2.py` is on `main` with per-phase `phase <name>: t … s, RSS …` prints and the build-only knob `FEM_EM_WF12_STEP4_BUILD_ONLY` (set only for the landing window; **it stays unset in this entry's command**, ruling 312-7). At `-n 2` the build-only window read (A0) 369 387 cells rel 0 and (Deg) 2 570 024 = 2·500 568 + 2·784 444 exactly, build 114.6 s, peak RSS 1.13 GiB summed (`20261008T093243Z_WF-12-step4a-W1.log`). The command is the review's to issue (below). *Was:* PENDING PREREQUISITE — §10 H16a (`WF-12` step 4a): the step-4 module `tests/validation/test_wf12_f_human_copper_degree2.py` landed on `main` from `attempt/WF-12-step4-20261006T1233Z` (`d3cd9bd8`) with per-phase time and RSS prints and a build-only knob, its (A0) 369 387 cells and (Deg) 2 570 024 = 2E + 2F asserted at `-n 2` as item 264 read them (`20261006T123338Z_WF-12-step4-W1.log`). **The implementer landing H16a marks this entry READY in the same commit and does not write the command**; the next review issues it in entry 30's durable-capture shape with that module, `FEM_EM_WF12_STEP4=1`, `-n 16`, timeout 14 400 s, and queues it for **Tuesday 2026-10-13 02:00** (the free night; the trailing-7-day count on 10-13 is six with it).

**Why:** the degree-2 column is the price of every tuned human-scale figure on copper (the full copper 32×32 through one factorisation, `PORT-19`; H17's root; H9′). Item 264 killed it at 590 s at `-n 8` with the factorising drive unreturned, so the ordinary tier cannot hold it and nothing else prices it.

**Price, measured in the family (scaling stated):** entry 21 (`POST-7` 2, the σ = 800 *solid* F-human at degree 2, 3 256 418 unknowns, `-n 16`) factorised in 303.8 s and ran one drive in 448.96 s at 104.39 GiB `memory.peak`; this hole mesh carries 2 570 024 unknowns (0.79×) but item 264 shows the copper wall + 32 sheets factor slower per unknown at `-n 8` (> 445 s for build + drive against the solid route's scaled ≈ 280 s). **Predicted 600–1 800 s, 60–110 GiB** at `-n 16`. Timeout 14 400 s; 512 GiB limit.

**Readout (record; `TH-14`'s copper identity and the exact power partition asserted at their imported bands on the column):** per-phase time (build, factorisation, drive) and RSS; the d1 → d2 move of the self and adjacent classes of column 17 against `WF-12` step 2's degree-1 record; the coil share.

**Decision rule (pre-registered; the 10-14 weekly reads it):** the identities green **and** the factorising drive ≤ 1 500 s ⇒ the full copper 32×32 at degree 2 (one factorisation + 31 held, `PORT-19`) is priced at ≈ (factorisation + 32 × held) and the 10-14 weekly commissions it as `xl` entry 34; > 1 500 s or memory > 200 GiB ⇒ the degree-2 copper F-human is an `xxl` case and is commissioned there or not at all; an identity red ⇒ a known-issues entry (the first degree-2 red on a copper hole mesh) and no 32×32. No band moves.

**Command** (issued 2026-10-08 10:30 review: entry 30's durable-capture shape with the module, `FEM_EM_WF12_STEP4=1` (`test_wf12_f_human_copper_degree2.py:45`, `:303`), the build-only knob unset (ruling 312-7), `-n 16`, and the raw-log name changed):

```
XL_CHUNK="WF-12-step4-H16"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/wf12-h16-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF12_STEP4=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf12_f_human_copper_degree2.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`.

### 33. `xl` — `IMP-3` step 7 / I15 part 3: the r/8 cap-ball rung at degree 2 on the exact-rectangle 1 g rod fixture, `-n 16` (weekly 2026-10-07, §10 chain I15)

**Status:** RUN `20261014T070008Z_IMP-3-step7-I15-r8.log` — 13 passed on all 16 ranks, Status 0, 511 s, 80.55 GiB `memory.peak`, 342 072 cells; I14's identities green at degree 2 on the r/8 rung; η_rod(T1) P1 3.482328 → 3.556262 → 3.530461 (Δ₁ +0.073934, Δ₂ −0.025801, signs opposite), printed (ledger row 2026-10-14, marked by the 2026-10-14 03:00 review; the decision rule below is the 10-14 weekly's). *Was:* QUEUED 2026-10-08 (10:30 review) for Wednesday 2026-10-14 02:00 — `docs/testing/xl-queue.d/20261014-IMP-3-step7-I15-r8.env` (command issued below by this review and copied verbatim; trailing-7-day count on 10-14 is six with it; the r/8 path's inputs `logs/imp3/ours_imp3g_d2off.json`, `ours_imp3g_d2rod.json`, `ours_imp3h_r4s2.json` are present and its output `ours_imp3h_r8.json` absent at queueing). *Was:* READY 2026-10-08 (§9 item 327, 07:30 slot). `FEM_EM_IMP3_STEP7` now accepts `r8` (the solve) and `r8build` (build only, landing window); on those two paths only `test_imp3_tip_ratio_1g.ROD_CELL_GUARD_W2` is monkeypatched to `ROD_CELL_GUARD_XL` = 400 000 (an `xl` resource guard, not a band), printed as a DISCLOSURE line and undone at teardown; (N2) checks the unpatched 320 000 before every dispatch. At `-n 2` the r8build window read (Kc8) `tip_size_m` 0.000125 == `TIP_R8` bitwise, 342 072 cells > 276 088; (G1x) 0.855 of the guard; (N1) 1.069 × 320 000; (G3) 8 790 tag-402 facets all exterior; (Rb8) "3.556262" string-equal; mesh identity with GEO-36 W2 exact (cells, 402 facets, 2 200 500 = 2E + 2F, near-T1 26 / 5.081837e-04 m, N_ball 93 132); build 83.57 s (`20261008T123452Z_IMP-3-item327-W1b-r8build.log`). `logs/imp3/ours_imp3h_r8.json` is absent (the r8 path writes it with mode `"x"` and checks absence before the build). **r/8 needs `xl` for the 590 s heavy ceiling, not for memory:** ≈ 67–72 GiB summed is predicted at `-n 8`, scaled from item 305 W1's measured 53.85 GiB at 276 088 cells. The n-port sweep honours `FEM_EM_SOLVER_PROGRESS=2` (read in code, not run: the first drive's `TimeHarmonicSolver.solve` sets `mat_mumps_icntl_4`, `time_harmonic.py:580–583`; drives 2–4 go through `solve_with_held_factorization`, which prints the back-substitution time, `:715–723`). The command is the review's to issue. *Was:* PENDING PREREQUISITE — **still pending after item 305 (2026-10-08 03:00 review):** item 305 landed the module green (`e6450a72`, W1 459 s, not killed, so its rung does not join this window), but its knob `FEM_EM_IMP3_STEP7` accepts only `r4s2` (module `:7–8`, `:243`) and I14's `_arm_d` asserts (G1) rod cells < `ROD_CELL_GUARD_W2` 320 000 against r/8's 342 072 — the rung is not yet selectable. **§9 item 327 (I15 part 3a) adds it; the implementer landing 327 marks this entry READY.** *Was:* PENDING PREREQUISITE — §9 item 305 (I15 part 2, the r/(4√2) rung, ruled heavy tonight) landing its module `tests/validation/test_imp3_tip_ladder_1g.py` on `main` with the rung selected by knob. **The implementer landing item 305 marks this entry READY in the same commit and does not write the command**; the next review issues it in entry 30's durable-capture shape with that module, the r/8 knob value, `-n 16`, timeout 14 400 s, and queues it for **Wednesday 2026-10-14 02:00** (the trailing-7-day count on 10-14 is six with it). If item 305 itself is killed at 590 s, its rung joins this window (two knob values, two pytest invocations in one command) and the review says so here.

**Why:** the tip ratio η_rod(T1) is not p-stable at 1 % on either mesh (−4.961 % record mesh, −5.995 % 1 g rung; I13 / I14) while the 1 g ball ratio moves 0.28 %; the question is discretisation at the tip, and the three-rung cap-ball ladder (r/4, r/(4√2), r/8) is the axis the quantity lives on (the 2026-09-30 epitaph's rule).

**Price, measured in the family (scaling stated):** I14's W3 (the r/4 rod arm at degree 2, 1 g rung, `-n 8`) ran 450 s at 50.97 GiB summed RSS (`20261006T172542Z_IMP-3-item284-W3-d2rod.log`); part 1 measured the cells per rung (252 118 → ≈ 290 k → 342 072) and scaled the rod arm to 603–621 s / 74–77 GiB at `-n 8`. **Predicted 400–900 s, 75–110 GiB at `-n 16`.** 342 072 cells exceed `ROD_CELL_GUARD_W2` (320 000): the window's module raises the guard for the knobbed rung only, disclosed, never for the default path. Timeout 14 400 s; 512 GiB limit.

**Readout (record; I14's identities (G1)–(G3), (i), (ii), (PP), (M), (E), (V), (D) asserted at their imported bands on the rung):** η_rod(T1) (P1, degree 2) at r/8 beside r/4's 3.482328 and item 305's r/(4√2) value; the measured near-T1 edge ratios per rung (≈ 1.087 then ≈ 1.264 in part 1 — unequal, so no observed order is claimed); the 1 g tip ratio Q at each rung.

**Decision rule (pre-registered; the 10-14 weekly reads it):** with Δ₁ = η(r/(4√2)) − η(r/4) and Δ₂ = η(r/8) − η(r/(4√2)): same sign and \|Δ₂\| ≤ 0.7 \|Δ₁\| ⇒ `IMP-3`'s row carries the three-rung bracket on η_rod(T1) as a record ("moving monotonically toward a limit; bracket [η(r/8) − \|Δ₂\|, η(r/8) + \|Δ₂\|]"), no rate; otherwise ⇒ the tip is not resolved by the cap ball alone (the near-T1 edge sits outside it) — a known-issues finding and §10 I17 (the ball-radius knob) opens. No band moves either way.

**Command** (issued 2026-10-08 10:30 review: entry 30's durable-capture shape with the module, `FEM_EM_IMP3_STEP7=r8` (`test_imp3_tip_ladder_1g.py:83–84`; not `r8build`), `-n 16`, and the raw-log name changed):

```
XL_CHUNK="IMP-3-step7-I15-r8"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/imp3-i15-r8-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_IMP3_STEP7=r8 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the r/8 path writes the gitignored `logs/imp3/ours_imp3h_r8.json`.

### 34. `xl` — `ANS-6` step 5c: the exact hole route at h = 0.005, degree 2, Cu and PEC, 64 MHz (weekly 2026-10-10, §10 chain F31)

**Status:** READ 2026-10-17 (weekly, §10 T41) — the copper-loss partition at 64 MHz converges: P_coil/P_in 0.4760 → 0.4611 → 0.4576 on h = 0.015 / 0.0075 / 0.005 (`:4894`; ratio 0.23 against first-order 0.33 / second-order 0.19 on those steps), limit ≈ 0.454, remaining ≤ 0.8 %; the Cu − PEC gap on the rung 0.0611 / 0.0277 / 0.0475 % against ×1's 0.0633 / 0.0287 / 0.0495 % (`:4906`). Read beside T30 as "ours alone": Q and `R_in` converge the same way ⇒ `ANS-7` ✅. *Was:* RUN `20261015T070008Z_ANS-6-step5c.log` — 25 passed, 14 skipped, Status 0, 6 928 s, 274.04 GiB `memory.peak`, 523 880 cells, 3 467 542 degree-2 unknowns (`:2425`, `:4898`); readout in the ledger row (2026-10-15 10:30 review); the reading is the 10-17 weekly's. *Was:* QUEUED 2026-10-11 (03:00 review) for Thursday 2026-10-15 02:00 — `docs/testing/xl-queue.d/20261015-ANS-6-step5c.env` (command below, copied verbatim). *Was:* READY 2026-10-10 (weekly) — the knobs are entry 31's
(`FEM_EM_EXACT_ANS6_RESOLUTION`, `_DEGREE`, `TH14_STEP2_FREQS_MHZ`; entry 31
runs them at 0.0075 on 10-12) and the command below is entry 31's with the
resolution and the raw-log name changed. The clerk queues it for **Thursday
2026-10-15 02:00** (trailing-7-day count on 10-15: 10-09, 11, 12, 13, 14, 15 =
six).

**Why:** `ANS-7`'s h = 0.005 rung (entry 28) moved Cu `R_in` and Q by −1.2 %
*away* from AED while `f₀` converged (T13, 2026-10-10): the Leontovich loss on
the copper hole mesh is not h-converged at 0.005. The copper-loss partition's
three-rung sequence (×1, 0.0075 from entry 31, 0.005 from this window) on
`ANS-6`'s fixture — one factorisation per column, no sweep — is the cheapest
rung that says whether the loss is converging, and it is the comparand of
`ANS-6`'s Finite-Conductivity column (private).

**Price, measured in the family (scaling stated):** entry 28 — the same
523 880-cell exact hole mesh at degree 2, 3 467 542 unknowns, five Cu
factorisations at 1 173–1 430 s each, 272.35 GiB `memory.peak`, mesh ≈ 110 s;
entry 14 (the ragged 0.005 hole mesh, PEC / Cu / P1) 4 110 s at 266.79 GiB. Two
columns (Cu and PEC) plus the identity ⇒ **predicted 2 800–4 200 s,
260–300 GiB** at `-n 16`. Timeout 14 400 s; 512 GiB limit.

**Readout (record; `TH-14`'s imported gates and the copper surface-loss
identity asserted on the rung):** Cu and PEC class means at degree 2 beside
F19's ×1 values and entry 31's 0.0075 rung; the Cu − PEC move per class on
each rung; the copper surface-loss partition P_surf / P_acc on each rung.
Private: the rung against the Finite-Conductivity and PEC columns.

**Decision rule (pre-registered; §10 T30 / F27, read by the 10-17 weekly):**
with Δ₁ = x(0.0075) − x(×1) and Δ₂ = x(0.005) − x(0.0075) on the copper-loss
partition: same sign and \|Δ₂\| ≤ ½ \|Δ₁\| ⇒ the loss is converging and
`ANS-6`'s row carries the bracket [x(0.005) − \|Δ₂\|, x(0.005) + \|Δ₂\|] as
the copper h-statement (the h steps are unequal, 0.015 → 0.0075 → 0.005, so no
rate is claimed); larger or opposite ⇒ the Leontovich loss on this fixture is
not h-resolved at 0.005 — a known-issues finding, and no finer copper rung is
bought before the strip-edge mechanism is named. No band moves.

**Command** (entry 31's shape with the resolution and the raw-log name
changed; `FEM_EM_EXACT_ANS6_PART` unset ⇒ both columns and the identity on
one mesh):

```
XL_CHUNK="ANS-6-step5c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/ans6-step5c-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 FEM_EM_EXACT_ANS6_RESOLUTION=0.005 TH14_STEP2_FREQS_MHZ=64 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the knobbed path writes only the gitignored `_r0.005` JSON,
never F19's records.

### 35. `xl` — `WF-10` step 4: the tuned two-port drive on the copper F-human at degree 2, κ off, 64 MHz, `-n 16` (weekly 2026-10-10, §10 chain H20)

**Status:** READ 2026-10-17 (weekly) — **second branch** (mean \|B₁⁺\| −23.4 %, SAR10g(c\*) −40.8 % at the degree-1 node, attributed to detuning: Im Y_in +1.346e-02 S, per phantom watt ≤ 1.6 %): the human-scale tuned figures are quoted at degree 2 **at the degree-2 root only** (entry 41); this window's values stay a labelled record of the detuned node. *Was:* **RUN `20261016T070008Z_WF-10-step4-H20.log`** (2026-10-16 02:00; Status 0, 18 passed / 1 skipped, 1 325 s, 86.81 GiB; every gate green; the d1 → d2 moves −23.4 % mean \|B₁⁺\| and −40.8 % SAR10g(c\*) fire the rule's **second branch** — ledger row; read by the 2026-10-16 03:00 review, the 10-17 weekly rules). **Caveat for the weekly (log-pathologist ruling):** the move is the detuning at the degree-1 node (Im Y_in(P17, (1,+1j)) +1.346e-02 S, `:6836`; P_phantom/P_src −39.9 %); per phantom watt \|B₁⁺\| moves +0.40 % and SAR10g(c\*) +1.58 %, so the order move of the *field* is inside 5 % and the move of the *tuned operating point* is not. *Was:* QUEUED 2026-10-13 (10:30 review) for **Friday 2026-10-16 02:00** — `docs/testing/xl-queue.d/20261016-WF-10-step4-H20.env` (the command below, copied verbatim; trailing-7-day count on 10-16: 10-11, 10-12, 10-13 run, 10-14, 10-15 queued, 10-16 = six). All three prerequisites MET; entry 38 is **not** folded in (the fold is the 10-14 weekly's call and its default is 10-19 on its own night). *Was (2026-10-13 03:00 review):* PENDING PREREQUISITE — (i) MET 2026-10-11 (below); (ii) **MET 2026-10-13** — entry 32's factorising drive 287.1 s ≤ 1 500 s and 89.31 GiB ≤ 200 GiB (`20261013T070008Z_WF-12-step4-H16.log:6710`, `:7810`); **(iii) added by this review: §9 item 414 (`WF-10` step 4b)** — `test_wf10_f_human_tuned_sar10g.py:164–184` reads the degree-1 `NODE_JSON` and asserts (S_rep) at `REUSE_REPRODUCTION_RTOL` (`:349`) against the degree-2 `s_node`, a certain red (entry 32's S17,17 d1 → d2 move 3.550 %, `:6715`); item 414 makes the module read `_node_json_path(_degree_knob())` and stop BLOCKED when the `-d2` record is absent. **MET 2026-10-13 (§9 item 414, the commit that adds this sentence; `20261013T110954Z_WF-10-step4b-W0.log`, `20261013T111006Z_WF-10-step4b-W1-regate.log`)**. **The implementer landing item 414 marks (iii) MET in the same commit and does not touch the command**; the next review then marks READY and queues it for **Friday 2026-10-16 02:00** (trailing-7-day count on 10-16: 10-11 … 10-16 = six). The command below is written by this review: entry 32's durable-capture shape, both step-2 modules in order, `FEM_EM_WF10_STEP2=1` added (without it the b1_sar gated test skips, `test_wf10_f_human_tuned_b1_sar.py:603–605`, and no `-d2` record is written — the 10-10 text named only `_STEP2_SAR`). **Price rewritten from entry 32:** each module calls the builder, so the window runs two builds and two degree-2 factorising node pairs (≈ 2 × (122 + 290) ≈ 830 s) plus 39 degree-2 10 g balls (unpriced at degree 2; ≈ 3.4 s each at degree 1) ⇒ **predicted 900–1 800 s, 90–190 GiB** (the upper figure if the first module's factorisation is not freed before the second build). Timeout 14 400 s; 512 GiB limit. *Was:* PENDING PREREQUISITE — (i) §10 H20a (`WF-10` step 4a): **MET 2026-10-11 (§9 item 368, the step 4a commit; `20261011T140313Z_WF-10-step4a-W1.log`, `20261011T140523Z_WF-10-step4a-W2-regate.log`)** — a degree
knob `FEM_EM_WF10_DEGREE` on the step 2 module
(`tests/validation/test_wf10_f_human_tuned_b1_sar.py` and its 10 g sibling),
build-only at `-n 2` asserting unknowns = 2E + 2F on the hole mesh and unset
bitwise on the degree-1 path; (ii) entry 32's reading (10-13): the factorising
degree-2 copper drive ≤ 1 500 s and memory ≤ 200 GiB. **The implementer
landing H20a marks (i) met in the same commit and does not write the
command**; the 10-14 03:00 review (after entry 32's log) marks READY and
issues the command in entry 32's durable-capture shape with the two step-2
modules, `FEM_EM_WF10_DEGREE=2`, `FEM_EM_WF10_STEP2_SAR=1`, `-n 16`, and
queues it for **Friday 2026-10-16 02:00** (trailing-7-day count on 10-16:
10-11 … 10-16 = six). If entry 32's drive exceeds 1 500 s, this entry is
re-tiered beside entry 36. At degree 2 the builder writes `logs/wf10/h9p-node-118.2437-d2.json` and 358's test runs `with_z50=False` ((R) / (Z50) are degree-1 comparands, printed only), so the command lists `test_wf10_f_human_tuned_b1_sar.py` before `test_wf10_f_human_tuned_sar10g.py`.

**Why:** H9′'s readings are degree 1 on a fixture whose degree-1 → 2 move on
F-small was 2–7 % on the self class; the human-scale tuned figures cannot be
quoted at degree 1 until their order move is measured once.

**Price (to be measured by entry 32; the family's scaling):** entry 32 prices
one factorising degree-2 copper drive at `-n 16` (predicted 600–1 800 s,
60–110 GiB); this window is one factorisation plus the two feed drives held
(`PORT-19`) and the `\|B₁⁺\|` / 10 g post-processing ⇒ predicted ≈ 1.2 × entry
32's drive, memory ≈ entry 32's. Timeout 14 400 s; 512 GiB limit. The
prediction is rewritten from entry 32's measured figures before READY.

**Readout (record; step 2's gate set asserted at degree 2 — the feed-bisector
mirror on `\|B₁⁺\|` and 10 g SAR at `C4_COVARIANCE_BAND`, (P_w), (P_RC),
(Σ_P), coverage, containment):** the d1 → d2 move of mean `\|B₁⁺\|_ccw`,
SAR10g(c\*), η and the k = ±1 share; Im `Y_in` at the degree-1 root (the
degree-2 root is not re-found here — its sign says which way it moved).

**Decision rule (pre-registered; the 10-17 weekly reads it):** every gate
green and the d1 → d2 moves of mean `\|B₁⁺\|` and SAR10g(c\*) ≤ 5 % ⇒
human-scale tuned figures may be quoted at degree 1 with the order move
beside them; a move > 5 % ⇒ they are quoted at degree 2 only and the
degree-2 root (H17 at degree 2, an `xl` family) is the next step; a gate red
⇒ a known-issues entry (the first degree-2 red on the tuned copper F-human).
No band moves.

**Command** (written 2026-10-13 03:00 review in entry 32's durable-capture shape; queued only after (iii) is MET):

```
XL_CHUNK="WF-10-step4-H20"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf10 && R=/workspace/logs/wf10-h20-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF10_DEGREE=2 FEM_EM_WF10_STEP2=1 FEM_EM_WF10_STEP2_SAR=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf10_f_human_tuned_b1_sar.py tests/validation/test_wf10_f_human_tuned_sar10g.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the degree-2 node record `logs/wf10/h9p-node-118.2437-d2.json` is gitignored.

### 36. `xxl` (or `xl`, tier decided by entry 32) — `WF-12` step 5: the degree-2 copper 32×32 on the F-human hole mesh, κ off, 64 MHz, `-n 16` (weekly 2026-10-10, §10 chain H16 → H)

**Status:** **READ 2026-10-19 (weekly, §10 rulings) — gates green ⇒ the degree-2 copper column set exists;** H25 read `C_mode`(d2, Cu) = 114.5221 pF from the saved matrix, −0.865 % of the root 115.5209 pF (the circuit layer's prediction on copper at degree 2 inside 1 % of the in-model root; its share / leak predictions missed, printed, recorded not re-posed); H30 reads the same matrix at the root after entry 41. *Was:* *(2026-10-18 10:30 review: the readout's "`ring_mode_rows_v2` on the saved matrix" clause is met by H25, §9 item 516, `20261018T095559Z_ANS-8-H25-W1.log` — 18 class means reproduced to ≤ 3.4e-10, `C_mode`(d2, Cu) 114.5221 pF.)* **RUN `20261018T070008Z_WF-12-step5-H16c.log`** (2026-10-18 02:00 `xl`; Status 0, 15 passed on all 16 ranks, 825 s, 87.93 GiB `memory.peak` — 2.3 % below the 90 GiB floor; every degree-2 gate green — (P) 9.4e-13, NEG ≥ 118.6×, (ii) 1.0e-13, (iii) σ_max 0.999983227, worst class spread 3.295 % of 5 %; the degree-2 32×32 saved to `logs/wf12/h16b-cu-s32-d2.npz`, bitwise reload — H25 is unblocked, §9 item 516; `xl-ledger.md` row, 2026-10-18 03:00 review). *Was:* QUEUED 2026-10-15 (03:00 review) for **Sunday 2026-10-18 02:00** — `docs/testing/xl-queue.d/20261018-WF-12-step5-H16c.env` (the command below, copied verbatim). *Was:* READY 2026-10-14 (weekly) — tier **`xl`** by the rule below (entry 32: factorising drive 287.1 s ≤ 1 500 s, 89.31 GiB ≤ 200 GiB), for **Sunday 2026-10-18 02:00** (trailing-7-day count on 10-18: 10-12, 13, 14, 15, 16, 18 = six); the command below was written tonight in entry 32's durable-capture shape with H16b's two knobs (`FEM_EM_WF12_STEP5=all`, the save path under the gitignored `/workspace/logs/wf12/`), `-n 16`, timeout 14 400 s; prerequisites (i) and (ii) both MET. The clerk copies it verbatim into `xl-queue.d/20261018-WF-12-step5-H16c.env`. *Was (2026-10-13 03:00 review):* PENDING PREREQUISITE — (i) entry 32's reading: **MET 2026-10-13** — factorising drive 287.1 s ≤ 1 500 s, 89.31 GiB ≤ 200 GiB ⇒ by the tier rule below this is an **`xl`** case (the 10-14 weekly re-tiers it); measured price ≈ 122 build + 287 factorising drive + 31 × 5.64 held ≈ 584 s at `-n 16` (`20261013T070008Z_WF-12-step4-H16.log:6762`), memory ≈ entry 32's 89 GiB plus the held columns' storage — **predicted 600–1 200 s, 90–120 GiB**; (ii) H16b — §9 item 415 (`WF-12` step 5a, the all-columns and save knobs, build-only at `-n 2`), opened by this review; the implementer landing it marks (ii) met and does not write the command (the 10-14 weekly writes it). `xxl` 10-17 is therefore not this entry's. *Was:* PENDING PREREQUISITE — entry 32's reading (10-13). **Tier
rule, pre-registered:** factorising drive ≤ 1 500 s and memory ≤ 200 GiB ⇒
the 10-14 weekly re-tiers this to `xl` (one factorisation + 31 held columns,
`PORT-19`; predicted ≈ factorisation + 31 × held + the identity); drive
> 1 500 s or memory 200–500 GiB ⇒ `xxl`, **Saturday 2026-10-17 02:00**;
memory > 500 GiB ⇒ not commissioned (the degree-2 copper F-human is
re-posed on the quarter, `TH-20` 2b). The 10-14 weekly writes the command in
entry 29's shape (`FEM_EM_WF12_STEP4=1` plus the all-columns and save knobs
the step-5 module adds — H16b, a lane-1 landing step the daily opens once
entry 32 reads; **H16b MET 2026-10-13 (§9 item 415, `a1b918b0`; `20261013T111749Z_WF-12-step5a-W1.log`)** — the knobs are `FEM_EM_WF12_STEP5=all` and `FEM_EM_WF12_STEP5_SAVE_NPZ=/workspace/logs/wf12/h16b-cu-s32-d2.npz` on `tests/validation/test_wf12_f_human_copper_s32_degree2.py`, which does not read `FEM_EM_WF12_STEP4`) and the measured price.

**Why:** the full copper 32×32 at degree 2 is the matrix every human-scale
tuned figure on copper is read from (H15c's degree-2 `C_mode` on copper,
`ANS-8` step 2's degree-2 column, H17 at degree 2); `xxl` 09-26 proved the
σ = 800 solid twin fits an `xl` night, and copper's price is entry 32's.

**Readout (record; `TH-14`'s copper identity, reciprocity, passivity and the
18 C16 × mirror classes asserted at their imported bands; the matrix saved
under the gitignored `/workspace/logs/wf12/`):** per-phase time and RSS; the
d1 → d2 move of the 18 class means against `WF-12` step 2's degree-1
record; `ring_mode_rows_v2` on the saved matrix (`C_mode`, k-share) beside
the degree-1 and σ = 800 rows.

**Decision rule (pre-registered):** gates green ⇒ the human-scale copper
degree-2 column set exists and `ANS-8` step 2's degree-2 column is read
from it; a class spread or passivity red at degree 2 ⇒ a known-issues entry
and the quarter route (`TH-20`) carries degree 2 instead.

**Command** (written 2026-10-14 weekly: entry 32's durable-capture shape with the step-5 module, `FEM_EM_WF12_STEP5=all` and `FEM_EM_WF12_STEP5_SAVE_NPZ` (`test_wf12_f_human_copper_s32_degree2.py:49–50`; the module does not read `FEM_EM_WF12_STEP4`), `-n 16`, and the raw-log name changed):

```
XL_CHUNK="WF-12-step5-H16c"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf12 && R=/workspace/logs/wf12-h16c-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF12_STEP5=all FEM_EM_WF12_STEP5_SAVE_NPZ=/workspace/logs/wf12/h16b-cu-s32-d2.npz timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf12_f_human_copper_s32_degree2.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the matrix is saved under the gitignored `/workspace/logs/wf12/`.

### 37. `xl` — `TH-20` step 2b / H13b: the four degree-2 quarter sectors reassembled against a degree-2 mirror-built P17 column solved in the same window (weekly 2026-10-12; re-worded 2026-10-14, §10 chain H13b)

**Status:** **RUN `20261020T070008Z_TH-20-step2b-H13b.log`** — 14 passed on all 16 ranks, Status 0, 936 s, 90.24 GiB (both inside the brackets); (A) 1.739e-11 against 1e-9 ⇒ the decision rule's first branch: `TH-20` ✅ (row written by the 2026-10-20 03:00 review; C1 / C1x / C2 printed, not asserted — the Readout's parenthetical, ruled non-gating, `TH-20.md`). *Was:* QUEUED 2026-10-15 (03:00 review) for **Tuesday 2026-10-20 02:00** — `docs/testing/xl-queue.d/20261020-TH-20-step2b-H13b.env` (the command below, copied verbatim). *Was:* READY 2026-10-14 (weekly) for **Tuesday 2026-10-20 02:00** (trailing-7-day count on 10-20: 10-14, 15, 16, 18, 19, 20 = six). Re-worded tonight as the 10-13 review asked: the comparand is **the degree-2 mirror-built P17 column solved in the same window** (not entry 32's direct-mesh column — the two meshes differ by 3e-4 even at degree 1), the module is step 1's `tests/validation/test_th20_quarter_solves.py` with `FEM_EM_TH20_STEP1=1 FEM_EM_TH20_STEP1_FULL_DRIVES=1 FEM_EM_TH20_DEGREE=2` (H13b-a landed the knob, §9 item 412); price rewritten from entry 32 and step 2a — four degree-2 sectors (684 904 unknowns each, 66–86 s per sector at `-n 8`, step 2a) plus the degree-2 mirror-built full column (2 702 016 unknowns, ≈ entry 32's 287 s drive at 1.05× the unknowns) ⇒ **predicted 500–1 000 s, 90–120 GiB** at `-n 16`; timeout 14 400 s; 512 GiB limit. The clerk copies the command below verbatim into `xl-queue.d/20261020-TH-20-step2b-H13b.env`. *Was (2026-10-13 03:00 review):* PENDING PREREQUISITE — entry 32 read: it did **not** leave the degree-2 column 17 on disk (it reads only the degree-1 `h12c-cu-s32.npz`, `20261013T070008Z_WF-12-step4-H16.log:6716`), so the withdrawal branch below does not fire and the entry stays for the 10-14 weekly. Prediction rewritten from entry 32's figures by this entry's own scaling: ≈ 0.65 × 287 s ≈ 190 s of sector solves plus the builds, memory ≈ 0.25 × 89 GiB ≈ 22 GiB. ~~That memory would fit the heavy tier, so the reading for the weekly is: once entry 36 runs with H16b's save knob (§9 item 415) the direct degree-2 column exists on disk and H13b can be a heavy lane-1 item.~~ *Retracted in the same review (item 412's drafter):* (A) compares the quarter against the **mirror-built** mesh (389 656 cells). Entry 32 / 36 solve the **direct** hole mesh (369 387 cells), and the two differ by 3e-4 on S(P17, P17) even at degree 1 (`test_th20_quarter_solves.py:60–62, :528`; `TH-20.md` ruling (3)). So no direct-mesh column can be (A)'s comparand. The degree-2 mirror-built column needs `xl`, which keeps this entry `xl`, and the withdrawal clause above cannot fire as written. **For the 10-14 weekly:** re-word the comparand ("entry 32's column 17" should read "a degree-2 mirror-built column solved in the same window"), and note the module: `tests/validation/test_th20_quarter_f_human.py` does not exist. The degree knob lands on step 1's `tests/validation/test_th20_quarter_solves.py` (§9 item 412, H13b-a, `FEM_EM_TH20_DEGREE`), so this entry becomes one window with `FEM_EM_TH20_STEP1_FULL_DRIVES=1 FEM_EM_TH20_DEGREE=2`. *Was:* PENDING PREREQUISITE — entry 32's reading (10-13). If entry 32 leaves the degree-2 column 17 on disk under the gitignored `/workspace/logs/wf12/` and the quarter solve at degree 2 fits the heavy ceiling (`TH-20` step 2a's price scaled: four sectors of ≈ 642 506 unknowns each, i.e. ≈ 0.16× entry 32's factorisation by N^1.33, ×4), H13b is a **lane-1 / heavy** item and this entry is withdrawn; otherwise the 10-14 weekly issues the command in entry 32's durable-capture shape with `tests/validation/test_th20_quarter_f_human.py`'s degree knob (H13b's landing step adds it, build-only at `-n 2`, as H20a did for `WF-10`) and the price measured from entry 32, for **Sunday 2026-10-18 02:00** (trailing-7-day count on 10-18: 10-12 … 10-16 and 10-18 = six). **H13b-a MET 2026-10-13 (§9 item 412, the commit that adds this sentence; `20261013T095500Z_TH-20-step2ba-W1.log`, `20261013T095737Z_TH-20-step2ba-W2-regate.log`)** — the knob is `FEM_EM_TH20_DEGREE` on `tests/validation/test_th20_quarter_solves.py` (no `test_th20_quarter_f_human.py` exists); with `FEM_EM_TH20_STEP1=1 FEM_EM_TH20_STEP1_FULL_DRIVES=1 FEM_EM_TH20_DEGREE=2` the module runs four degree-2 sectors and the degree-2 mirror-built P17 column and asserts (A) on its 32 entries.

*2026-10-21 weekly:* **READ** — the 10-20 03:00 review's close-out ratified: (A) 1.739e-11 against 1e-9 on 14 passed at 16 ranks, 936 s / 90.24 GiB inside both brackets ⇒ `TH-20` ✅, the quarter the licensed degree-2 human-scale route. The printed C1 / C1x / C2 controls get no control window (no consumer). No re-run.

**Why:** `TH-20`'s Done-when (ruled 2026-10-10) needs (A), the degree-2 reassembly of the mirror-built F-human copper 32×32's column against a direct degree-2 column; step 1 did it at degree 1 to 2.1e-11. The quarter is also the route that carries degree 2 at human scale if entry 36's memory exceeds 500 GiB.

**Price (to be measured by entry 32; the family's scaling):** step 1 ran four degree-1 quarter sectors plus the reassembly in 369 s (`20261005T020907Z_TH-20-step1-W2.log`); at degree 2 each sector carries ≈ 2 570 024 / 4 unknowns. Predicted ≈ 0.65 × entry 32's drive time in total, memory ≈ 0.25 × entry 32's peak. Timeout 14 400 s; 512 GiB limit. Rewritten from entry 32's figures before READY.

**Readout (record; step 1's identities C1 / C1x / C2 asserted at their imported bands at degree 2):** the reassembly residual (A) against the degree-2 mirror-built P17 column solved in the same window (band `1e-9`, step 1's), per-sector time and RSS, the full column's time and RSS.

**Decision rule (pre-registered):** (A) ≤ 1e-9 ⇒ `TH-20` ✅ and the quarter is licensed as the degree-2 human-scale route; a residual above it with the degree-1 reassembly still green ⇒ a known-issues entry (the first order-dependent symmetry red) and `TH-20` stays 🟡. No band moves.

**Command** (written 2026-10-14 weekly: entry 32's durable-capture shape with step 1's module and the three knobs (`test_th20_quarter_solves.py:49–50`, `:91`), the build-only knob unset, `-n 16`, and the raw-log name changed):

```
XL_CHUNK="TH-20-step2b-H13b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/th20 && R=/workspace/logs/th20-h13b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_TH20_STEP1=1 FEM_EM_TH20_STEP1_FULL_DRIVES=1 FEM_EM_TH20_DEGREE=2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_th20_quarter_solves.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the saved R goes under the gitignored `/workspace/logs/th20/`.

### 38. `xl` — `MAT-15` step 6 / P22: the tuned copper F-human's heating at degree 2 (weekly 2026-10-12, §10 chain P22)

**Status:** **SUPERSEDED 2026-10-17 (weekly) by entry 42** — this command heats the phantom from the degree-1 node, which entry 35 showed is detuned at degree 2 (the module computes its own Z_RC from `C_IN_F_HUMAN_PF`, `test_mat15_f_human_tuned_heating.py:360`), so it would measure the matching term. The clerk replaces `20261019-MAT-15-step6-P22.env` with entry 42's command once §10 P22b lands; if P22b is not on `main` by the 10-18 18:00 review, the clerk swaps the 10-19 and 10-20 queue files (entry 37 runs Monday, entry 42 Tuesday). *Was:* QUEUED 2026-10-15 (03:00 review) for **Monday 2026-10-19 02:00** — `docs/testing/xl-queue.d/20261019-MAT-15-step6-P22.env` (the command below, copied verbatim). *Was:* READY 2026-10-14 (weekly) for **Monday 2026-10-19 02:00** on its own night (trailing-7-day count on 10-19: 10-13, 14, 15, 16, 18, 19 = six) — **not folded into entry 35** (35 is queued and its file is not edited; one factorisation per window is the simpler record). Runs after entry 35, so the degree-2 node record it may read exists. Prerequisites (i) P22a and (ii) entry 32's reading both MET (below); the command was written tonight in entry 35's shape with the step-5 module, `FEM_EM_MAT15_STEP5=1`, `FEM_EM_MAT15_STEP5_DEGREE=2` together with `FEM_EM_WF10_DEGREE=2` (the module raises otherwise), `FEM_EM_MAT15_STEP5_SERIES=1`, `-n 16`, timeout 14 400 s. The clerk copies it verbatim into `xl-queue.d/20261019-MAT-15-step6-P22.env`. *Was (2026-10-13 03:00 review):* PENDING PREREQUISITE — (ii) **MET** (entry 32: drive 287.1 s, 89.31 GiB); (i) P22a queued as §9 item 413. Entry 35's command, written by this review, carries the two `WF-10` modules only and waits on its own prerequisite (iii) (§9 item 414); if P22a lands before entry 35 is marked READY the fold below is still open, and the 10-14 weekly says which — by default this entry runs on its own night (10-19). Prediction rewritten from entry 32's figures: ≈ 1.2 × 287 + 122 build + 200 s CN ≈ 670 s, memory ≈ 90–100 GiB ⇒ **predicted 600–1 200 s, 90–110 GiB**. *Was:* PENDING PREREQUISITE — (i) §10 P22a (lane 1, no solve): **MET 2026-10-13 (§9 item 413, the commit that adds this sentence; `20261013T110156Z_MAT-15-step6a-W1.log`, `20261013T110401Z_MAT-15-step6a-W2-regate.log`)** — `FEM_EM_MAT15_STEP5_DEGREE=2` must be set together with `FEM_EM_WF10_DEGREE=2` (the module raises otherwise); at degree 2 the series goes to `logs/mat15_s5/theta-series-ccw-d2.json` and (S_rep) is printed only — a degree knob on `MAT-15` step 5's module (`FEM_EM_MAT15_STEP5_DEGREE`, H20a's pattern — build-only at `-n 2` asserting unknowns = 2E + 2F on the hole mesh, unset bitwise on the degree-1 path); (ii) entry 32's reading (10-13): the factorising degree-2 copper drive ≤ 1 500 s and memory ≤ 200 GiB. The implementer landing P22a marks (i) met in the same commit and does not write the command; the 10-14 weekly issues it in entry 35's shape with the step-5 module, the degree knob and `FEM_EM_MAT15_STEP5_SERIES` set, `-n 16`, for **Monday 2026-10-19 02:00** (trailing-7-day count on 10-19: 10-13 … 10-16, 10-18, 10-19 = six) — or folds it into entry 35's command as a third module if P22a lands before entry 35 is marked READY (one factorisation serves both; the 10-14 weekly says which).

**Why:** P17 is the only thermal figure on a tuned coil at human scale and it is degree 1 on the 0.015 m phantom; `MAT-20` step 3's dose and six-minute figures are read from its series. The d1 → d2 move of the ball-mean θ(t) and of the 10 g SAR source is the record that licenses (or not) quoting them at degree 1, exactly as entry 35 does for `\|B₁⁺\|`.

**Price (to be measured by entry 32; the family's scaling):** entry 35's (one factorisation plus the two feed drives held) plus the CN thermal solve, which is P1 on the phantom and cost 160 s at `-n 8` in item 380's window. Predicted ≈ 1.2 × entry 32's drive + 200 s, memory ≈ entry 32's. Timeout 14 400 s; 512 GiB limit.

**Readout (record; `MAT-15` step 5's identities asserted at their imported bands — the energy identity, the calorimetric slope, the mirror on the 10 g source):** the d1 → d2 move of ⟨θ⟩(t\*), the 10 g ball-mean θ(t\*) at c\*_θ and SAR10g(c\*) per accepted watt; `cem43` and the six-minute maximum on the degree-2 series beside item 380's.

**Decision rule (pre-registered; the 10-21 weekly reads it — the 10-17 weekly precedes the window):** every identity green and the d1 → d2 move of the ball-mean θ(t\*) ≤ 5 % ⇒ P17's figures may be quoted at degree 1 with the order move beside them; > 5 % ⇒ degree 2 only; an identity red ⇒ a known-issues entry. No band moves.

**Command** (written 2026-10-14 weekly: entry 35's durable-capture shape with the step-5 module and its knobs (`test_mat15_f_human_tuned_heating.py:149–154`; the degree knob is guarded by `FEM_EM_WF10_DEGREE`, `:171–180`), `-n 16`, and the raw-log name changed):

```
XL_CHUNK="MAT-15-step6-P22"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf10 /workspace/logs/mat15_s5 && R=/workspace/logs/mat15-p22-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF10_DEGREE=2 FEM_EM_MAT15_STEP5=1 FEM_EM_MAT15_STEP5_DEGREE=2 FEM_EM_MAT15_STEP5_SERIES=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_mat15_f_human_tuned_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the degree-2 series goes to the gitignored `logs/mat15_s5/theta-series-ccw-d2.json`.

### 39. `xxl` (Saturday 2026-10-17, if READY by the 10-16 03:00 review) or `xl` (Wednesday 2026-10-21) — `WF-10` step 1e at degree 2 / H26: the in-model Larmor root on the copper F-human with r = 0.2 Ω at degree 2, `-n 16` (weekly 2026-10-14, §10 chain H26)

**Status:** READ 2026-10-17 (weekly, §10 rulings) — **first branch fires:** one − → + crossing on all four readings in [115.5, 116.0] pF, C_in(d2) = 115.520946 pF (four-reading spread 3.1e-4 rel), the grid-miss flag does not arise (the 118.2437 node re-measured entry 35's sign and value); H9′'s readings are re-posed at the root — entry 41 READY, re-tiered to `xl` Thursday 2026-10-22; `ans:8` v3 rows 3–5 at degree 2 are §10 H30, serial on entry 41. *Was:* **RUN `20261017T070008Z_WF-10-step1e-H26.log`** (2026-10-17 02:00 `xxl`; Status 0, 17 passed / 1 skipped on all 16 ranks, 4 018 s, 91.09 GiB `memory.peak`; every degree-2 gate green; one − → + crossing on all four readings in [115.5000, 116.0000] pF, C_in(d2) = 115.520946 pF — the first branch reads as firing, the weekly's to rule; the root is a secant estimate with no node record, and the downstream modules read a fixed 118.2437 pF node, so H27 at the root needs a node knob first; `xxl-ledger.md` row, 10:30 review). *Was:* **QUEUED 2026-10-16 (03:00 review) for Saturday 2026-10-17 02:00 (`xxl`)** — `docs/testing/xxl-queue.d/20261017-WF-10-step1e-H26.env`, the command under "Command as queued" below (`xxl` trailing-7-day count on 10-17: none since 10-03 — one of one). **READY 2026-10-16 03:00:** (i) MET 10-15 (below); (ii) **MET** — entry 35 ran 10-16 (`20261016T070008Z_WF-10-step4-H20.log`, Status 0) and its d1 → d2 moves (mean \|B₁⁺\| −23.4 %, SAR10g(c\*) −40.8 %) fire the second branch, confirmed by a `log-pathologist` ruling, which adds that the move is the detuning at the degree-1 node (per phantom watt ≤ 1.6 %) — the case *for* finding the degree-2 root. **The 15:00 notes (a)–(c) are settled here:** (a) the grid — Im Y_in(P17, (1,+1j)) at 118.2437 pF is **+1.346e-02 S** at degree 2 (`:6836`) and Im Y_in rises with C at degree 1 (every interval slope positive, `20261010T094215Z_WF-10-step1e2-W2.log:65`), so the root lies **below** 118.2437 pF; the degree-2 reading sits beside the degree-1 reading at 121.0000 pF (`:64`), a ≈ −2.75 pF translation, so the root is extrapolated near 115.5 pF, ≈ 115.4–116.3 pF over the degree-1 slopes, below 114 pF only if the degree-2 slope falls under ≈ 3.2e-3 S/pF. The ruling warns that this is a one-point extrapolation, so the fixed eight-node grid is **112.5000, 114.0000, 115.0000, 115.5000, 116.0000, 116.5000, 117.2500, 118.2437 pF**. It reaches 2.9 pF below the extrapolation and keeps the measured-positive 118.2437 node as its upper end (that node re-measures entry 35's reading in this module, a printed cross-check). Its resolution is 0.5 pF (≈ 4.3e-3 rel) across 115.0–116.5. (b) the Readout is read against the gates the module asserts at degree 2: (A0), (iv), (Deg), node kinds `["solve", "held"]`, (P50), NEG, (P_RC) per node (`test_wf10_f_human_copper_esr_root.py:596–607`), with (R) / (Z50) printed; (c) the Z50 pair is **kept** — it carries (P50) / NEG. **Price rewritten from entry 35:** build 145–161 s, Z50 pair ≈ 365 s, eight node pairs at 328–361 s ⇒ **predicted 2 800–4 500 s, 87–110 GiB** (entry 35 peak 86.81 GiB); timeout 14 400 s inside `xxl`'s 8 h. *Was:* PENDING PREREQUISITE — (i) §10 H26a (lane 1, no solve): **MET 2026-10-15 (§9 item 459, this entry's landing commit — `git log -1 -- tests/validation/test_wf10_f_human_copper_esr_root.py`; `20261015T093246Z_WF-10-H26a-W1.log`, `20261015T093452Z_WF-10-H26a-W2-regate.log`, `20261015T093834Z_WF-10-H26a-W3-regate.log`)** — a degree knob on `tests/validation/test_wf10_f_human_copper_esr_root.py` in H20a's pattern (`FEM_EM_WF10_DEGREE` read at call time; build-only at `-n 2` asserting unknowns = 2E + 2F = 2 570 024 on the hole mesh; the degree-1 path bitwise unchanged when unset); (ii) entry 35's reading (10-16 02:00): its rule's **second** branch — the d1 → d2 move of mean `\|B₁⁺\|` or SAR10g(c\*) **> 5 %** — fires; a move ≤ 5 % **withdraws this entry** (the degree-1 root stands with the order move beside it). **The implementer landing H26a marks (i) MET in the same commit and does not write the command**; the 10-16 03:00 review, after entry 35's log, marks READY if (ii) fires, checks the knob names below against H26a's landing commit, and queues it into `xxl-queue.d/20261017-WF-10-step1e-H26.env` for Saturday 02:00 (an `xl`-sized job on the `xxl` night is allowed; the tier is a ceiling). If (i) is not landed by then and (ii) fires, the clerk queues it for **Wednesday 2026-10-21 02:00** (`xl`; trailing-7-day count on 10-21: 10-15, 16, 18, 19, 20, 21 = six) once H26a lands. Knobs as landed: `FEM_EM_WF10_DEGREE=2` (H20a's parser, imported), `FEM_EM_WF10_STEP1E=1`, `FEM_EM_WF10_STEP1E_ROOT=1`; `FEM_EM_WF10_STEP1E_NODES_PF` is **required** (the root reader asserts it non-empty, `test_wf10_f_human_copper_esr_root.py:683`, and unset it gives one node); `FEM_EM_WF10_STEP1E_BUILD` stays unset; at degree 2 the nodes go to `logs/wf10/h17-node-<pF>-d2.json` and (R) / (Z50) are printed only.

**Why:** entry 35's rule says a degree 1 → 2 move above 5 % on the tuned figures quotes them at degree 2 only, and then the root they are read at must be the degree-2 root — H17's secant on C_in at degree 2 — or the figures are read at a node that is not the resonance at that order (`ANS-8`'s 10-03 INCONCLUSIVE was exactly this, 2.18 % of `f₀`).

**Price, measured in the family (scaling stated):** entry 32 — one degree-2 copper factorising drive 287.1 s at 89.31 GiB `memory.peak` (build 122 s); H17 at degree 1 found the root in two heavy windows with 5–7 secant nodes. One factorisation per node, memory freed between nodes ⇒ **predicted 1 500–3 000 s (6–8 nodes), 90–110 GiB** at `-n 16`. Timeout 14 400 s; 754 GiB (`xxl`) or 512 GiB (`xl`) limit.

**Readout (record; H17's identities — the (1, −1j) root's circuit ↔ field closure at `REDUCTION_BAND`, (P_w), the mirror — asserted at degree 2):** the degree-2 node C_in,F-human(d2) beside 118.2437 pF (d1), the sign and size of Im Y_in at the degree-1 node, the number of nodes, per-node time and RSS; the degree-2 node JSON saved under `/workspace/logs/wf10/`.

**Decision rule (pre-registered; the weekly after the window reads it):** one − → + crossing of Im Y_in inside the secant's own resolution ⇒ the degree-2 root exists and H9′'s readings are re-posed at it (H27 at the degree-2 node, `ans:8` v3 rows 3–5 at degree 2); no crossing in 8 nodes ⇒ a known-issues entry (the first degree-2 root failure on copper) and the tuned human-scale figures stay degree-1 records with the order move beside them. No band moves.

**Command** (provisional shape — entry 35's with the root module and its knobs `FEM_EM_WF10_STEP1E=1 FEM_EM_WF10_STEP1E_ROOT=1` (`test_wf10_f_human_copper_esr_root.py:47`, `:547`) plus `FEM_EM_WF10_DEGREE=2` from H26a; **the review that marks READY confirms every knob name against H26a's landing commit before copying**):

```
XL_CHUNK="WF-10-step1e-H26"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf10 && R=/workspace/logs/wf10-h26-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF10_DEGREE=2 FEM_EM_WF10_STEP1E=1 FEM_EM_WF10_STEP1E_ROOT=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf10_f_human_copper_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the degree-2 node JSON is gitignored.

**Command as queued** (2026-10-16 03:00 review). It is the provisional shape above with three changes, all of them the READY review's to make: the `xxl` profile and service, because it runs on the Saturday `xxl` night; `FEM_EM_WF10_STEP1E_NODES_PF` set to grid (a); and the knob names, checked against H26a's landing (`test_wf10_f_human_copper_esr_root.py:47–48`, `:617`; `FEM_EM_WF10_DEGREE` through `_degree_knob`, `:367–369`). Nothing else is changed:

```
XL_CHUNK="WF-10-step1e-H26"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf10 && R=/workspace/logs/wf10-h26-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF10_DEGREE=2 FEM_EM_WF10_STEP1E=1 FEM_EM_WF10_STEP1E_ROOT=1 FEM_EM_WF10_STEP1E_NODES_PF=112.5000,114.0000,115.0000,115.5000,116.0000,116.5000,117.2500,118.2437 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf10_f_human_copper_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

**Flag for the weekly that reads it (the pre-registered rule is unchanged; a daily does not re-word it):** with eight same-sign nodes, the rule's "no crossing in 8 nodes ⇒ a known-issues entry" cannot tell a root failure from a grid miss. If all eight read positive, the root sits below 112.5 pF (a C move of more than 4.9 %). If the 118.2437 node reads negative, it contradicts entry 35's +1.346e-02 S on the same fixture. The `log-pathologist` raised this risk before the window. The weekly decides which reading it is.

**For the 10-16 03:00 review, before any READY (2026-10-15 03:00 review, from H26a's brief drafter, §9 item 459; not a command edit):** (a) the provisional command sets no node list — without `FEM_EM_WF10_STEP1E_NODES_PF` the module solves one node and the root reader's assert (`test_wf10_f_human_copper_esr_root.py:612–613`) is a certain red; the module has no adaptive secant, so the READY review chooses a fixed 6–8-node grid (≤ 4 decimals in pF) on the side of 118.2437 pF that entry 35's sign of Im Y_in names; (b) the Readout names gates the root module does not assert (the (1, −1j) circuit ↔ field closure, (P_w), the mirror) — after H26a it asserts (A0), (iv), (Deg), kinds, (P50), NEG and (P_RC) at degree 2, so the Readout is read against those; (c) keeping the Z50 pair at degree 2 prices the window at ≈ 2 200 s (6 nodes) – 2 800 s (8), inside the 1 500–3 000 s bracket; dropping it saves ≈ 290 s and loses (P50) / NEG. If (a)–(c) cannot be settled by that review, entry 39 is not spent on `xxl` 10-17 and returns to the 10-17 weekly.

### 40. `xl` — `GEO-38` step 1b / C1 (b): one degree-1 copper column on F-clinical-1.5T r1, in-memory build, `-n 8` (daily licence: **cost probe**, 2026-10-16 03:00 review; §10 chain C1)

**Status:** **RUN `20261021T070008Z_GEO-38-step1b-C1b.log`** (2026-10-21 02:00; Status 0, 598 s, 40.37 GiB; ledger row filled by the 2026-10-21 03:00 review). **Read by that review (the first after the window):** every assertion green, so the decision rule's first branch fires — `GEO-38` step 1's Done-when is met and the row is ✅ as a probe (audit PASS). Brackets: t_b 317.0 s inside; t_col 199.0 s, t_held 2.03 s and the 598 s window **below** their brackets; `memory.peak` 40.37 GiB inside. C3 = 644 s by the module's projection ⇒ `xl` (entry 44's tier and command are the 10-21 weekly's, then the 10-22 03:00 review's as written). The printed degree-2 projection carries the F-human exponent, not the measured 1.433 — the weekly prices C5 from the measurement. *Was:* RATIFIED 2026-10-17 (weekly, §10 C1) as the daily wrote it; the 10-22 03:00 review reads the brackets, and C3 (entry 44) is `xl` or larger on any outcome. *Was:* **QUEUED 2026-10-16 (10:30 review) for Wednesday 2026-10-21 02:00** — `docs/testing/xl-queue.d/20261021-GEO-38-step1b-C1b.env` (the command below, copied verbatim; `cmp` identical). *Was:* READY (2026-10-16 04:30 slot, §9 item 494): `tests/validation/test_geo38_f_clinical_column.py` lifted onto `main` by path from `attempt/GEO-38-20261016T004700Z` (`1453764e`); W1 0 passed / 2 skipped (`20261016T093031Z_GEO-38-c1b-lift-W1-census.log`), W2 green, printed `64 1147366 1 64000000.0` (`20261016T093044Z_GEO-38-c1b-lift-W2-imports.log`). Command untouched.
*2026-10-21 weekly:* **READ** — the 03:00 review's reading ratified (`GEO-38` ✅ as a probe, audit PASS): every identity green; t_b 317.0 s, t_col 199.0 s, t_held 2.03 s, 598 s, 40.37 GiB at `-n 8`, degree-1 exponent 1.433. Read into entries 44 (C3 on the `xxl` night, `-n 16`; price from this column) and 45 (I34, brackets as pre-registered with the tighter prediction beside them), and into C5's degree-2 price (≈ 1 774 s at the measured exponent, `xxl`-only, after C4). The "> 228 s" t_col floor was a loaded-box kill time, not a bound. No re-run.
- **Intended run: Wednesday 2026-10-21 02:00** (`xl`). The trailing-7-day count on 10-21: 10-15 and 10-16 have run, 10-18, 10-19 and 10-20 are queued, and 10-21 makes six.
- Entry 39 went to `xxl` 10-17, so its `xl` fallback on 10-21 lapses.
- **The 10-17 weekly owns C1 (b)'s commission (§10 C1) and may ratify or withdraw this entry before the night.**

**Why:** C1 (b)'s target is the measured price of one degree-1 factorising copper column. Both heavy attempts died at 590 s inside the mesh build:
- item 469 built in 361.6 s at `-n 8`, so the column ran > 228 s and was never read (`20261015T172143Z_GEO-38-s1b-W1.log`);
- item 485's gmsh alone took 431 s at `-n 4` (`20261016T003426Z_GEO-38-s1b-a2-WA-cache.log`).

C3's tier and C5's degree-2 price are read from this column. It is unpriced at every tier, so its first window is a cost probe.

**Price (predicted; the readout is the price):**
- build t_b 361.6 s (item 469, `-n 8`);
- column 39.2·(U₁/500 568)² ≈ 364 s with U₁ = 1 525 193 (item 469's lower bound already rules out the ρ^0.84 reading of 102 s);
- held drive ≈ 13 s;
- post-processing and startup ≈ 92 s.

The window is ≈ 830 s. Summed RSS is 46.2 GiB (10.20·ρ^(4/3), ρ = 3.1061), or 55.8 GiB if it scales as ρ^1.5. Timeout 14 400 s; 512 GiB limit.

**Readout (record; fixed brackets, no gate, no band):**
- **Brackets:** t_b [300, 650] s; t_col [228, 1 200] s; t_held [5, 60] s; window [600, 2 400] s; `memory.peak` [30, 120] GiB. A reading outside its bracket is reported as such.
- **Printed:**
  - RSS per rank and summed after the column and at the end;
  - the measured exponent ln(t_col/39.2)/ln ρ;
  - MUMPS statistics (`FEM_EM_SOLVER_PROGRESS=2`, verbosity only);
  - the PROJECTIONS line.
- **Asserted, by item 469's module as written; this entry adds none:**
  - before the column: (A0) cells 1 147 366 within `CELL_COUNT_BAND`, the (iv) and outer/cavity census, and (S) cavity max < R_s;
  - after the column: (P) on P33 / P34 ≤ 1e-6, P_surf > 0, NEG, and (R1) ≤ `RECIPROCITY_BAND`.

PRICE prints before the post-column asserts, so a red there still records the price. Expected census: 5 passed, 1 skipped (`=cache` skips). The module does not assert the rank count; `-n 8` keeps item 469's.

**Decision rule (pre-registered; the next weekly reads it):**
- **Column measured and every assertion green:** `GEO-38` step 1's Done-when is met (census 64 by the 2026-10-15 ruling), and the review audits 🟡 → ✅ as a probe. C3's tier follows from C3 = t_b + t_col + 63·t_held under §10's "heavy if ≤ 590 s, else `xl`". Since t_b + t_col > 590 s on any outcome, this reads `xl`, and the weekly commissions C3 at this price.
- **Column measured with (P), (R1) or NEG red and (iv) green:** the price stands as a record, `GEO-38` stays 🟡, and the review opens a known-issues entry for the clinical-scale Leontovich wall or shield.
- **(A0), (iv) or (S) red:** no price; the r1 mesh differs on the `xl` image. The step returns to the weekly, with no second window bought.
- **Killed or OOM:** the last flushed phase line is the lower bound, and the step returns to the weekly.

**Command** (entry 35's durable-capture shape, route `=1` — one window, in-memory build plus column — at `-n 8`):

```
XL_CHUNK="GEO-38-step1b-C1b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/geo38-c1b-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_STEP1B=1 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_column.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the `=1` route writes no file.

### 41. `xxl` (Saturday 2026-10-24) or `xl` (Thursday 2026-10-22) — `WF-10` step 4 at the degree-2 root / H28′: the tuned two-port drive on the copper F-human at degree 2, κ off, 64 MHz, C_in(d2) = 115.5209 pF, `-n 16` (daily licence: **priced-family variant** of entry 35, 2026-10-17 10:30 review; §10 chain H28)

**Status:** **READ 2026-10-24 (weekly, §10 rulings) — second branch, as the 03:00 review read it: every gate green at the degree-2 root, the root-to-root moves −3.17 % / −6.30 % (SAR10g > 5 %) ⇒ the human-scale tuned figures are quoted at degree 2 only, at C_in(d2) = 115.5209 pF, from this window, with the per-phantom-watt reading (+0.60 / +1.13 %) printed beside them, never substituted (the rule was root to root); H31's carry ratified 37 / 37. Phase 6 EXITED on this reading (§6, §2.1). Price OUTSIDE-below on time (a record), INSIDE on memory.** *Was:* **RUN `20261022T070008Z_WF-10-step4-H28.log`** (2026-10-22 02:00; Status 0, 1 075 s, 87.12 GiB; 20 passed / 2 skipped on 16 ranks; every Readout gate green; \|Im Y_in(P17, (1,+1j))\| 1.0009e-05 S — at the root). The rule, read literally by the 10-22 03:00 review (ledger row): **second branch** — mean \|B₁⁺\| −3.17 %, SAR10g(c\*) **−6.30 %** root to root ⇒ quoted at degree 2 only, at C_in(d2); per phantom watt the moves are +0.60 % / +1.13 % (P_phantom/P_src −4.73 %). The decision is the 10-24 weekly's. *Was:* **QUEUED 2026-10-18 (03:00 review) for Thursday 2026-10-22 02:00** — `docs/testing/xl-queue.d/20261022-WF-10-step4-H28.env` (the command below, copied verbatim; the knob `FEM_EM_WF10_NODE=d2root` confirmed against item 511's landing `12df58ce`, `test_wf10_f_human_tuned_b1_sar.py` `NODE_ENV`). *Was:* **READY 2026-10-17 (weekly) — re-tiered to `xl` Thursday 2026-10-22 02:00** (trailing-7-day count on 10-22: 10-16 run, 10-18 / 19 / 20 / 21 queued, 10-22 makes six; the `xxl` 10-24 night goes to entry 44). (ii) **MET** by the weekly's reading of entry 39 (first branch; the re-pose uses C_in(d2) = 115.5209 pF, the landed `C_IN_F_HUMAN_D2_PF`). The clerk copies the command below (profile and service `xl`, as written) into `xl-queue.d/20261022-WF-10-step4-H28.env` after confirming the knob name against item 511's landing commit. *Was:* PENDING PREREQUISITE.
- **(i) MET 2026-10-17 (§9 item 511, the item-511 outcome commit on `main`; `20261017T170238Z_WF-10-H28-W1.log`, `20261017T170443Z_WF-10-H28-W2-regate.log`, `20261017T170806Z_WF-10-H28-W3-regate.log`).** Landed knob names: `FEM_EM_WF10_NODE=d2root` (with `FEM_EM_WF10_DEGREE=2`); constant `C_IN_F_HUMAN_D2_PF = 115.5209`; build-only selector `FEM_EM_WF10_STEP4C_BUILD` (W1 only, never in this entry). (ii) not yet MET, so the entry stays PENDING. §9 item 511 (H28): the node knob `FEM_EM_WF10_NODE=d2root` on `tests/validation/test_wf10_f_human_tuned_b1_sar.py` / `test_wf10_f_human_tuned_sar10g.py`. The implementer landing it marks (i) MET in the same commit, and also marks this entry READY if (ii) already reads MET. It never touches the command.
- **(ii) the 2026-10-17 weekly's reading of entry 39:** that its one − → + crossing at 115.520946 pF (four-reading spread 0.035723 pF, `20261017T070008Z_WF-10-step1e-H26.log:9000`) fires entry 39's first branch, "H9′'s readings are re-posed at it" (`:1877` of this file), and which node value the re-pose uses. A weekly that adopts entry 35's per-phantom-watt reading instead, with no re-pose, **withdraws this entry**.
- **Licence:** step 6b rule 4. The `xxl` floor is short (0 ahead after entry 39 ran); `xl` has 4 ahead (36, 38, 37, 40). The script, mesh, degree and rank count are entry 35's, and one named env knob is varied.
- **Night:** `xxl` Saturday 2026-10-24 02:00 by default. On 10-24 the `xxl` trailing 7 days hold nothing since 10-17, so this is one of one. An `xl`-sized job on the `xxl` night is allowed; the tier is a ceiling (entry 39's precedent). The weekly may re-tier it to `xl` Thursday 2026-10-22 02:00: on 10-22 the trailing 7 days hold 10-16 run and 10-18 / 19 / 20 / 21 queued, so 10-22 makes six.

**Why:** entry 35 measured the d1 → d2 move at the **degree-1** node. Most of that move was read as detuning (`log-pathologist`, 2026-10-16 03:00 review): Im Y_in(P17, (1,+1j)) was +1.346e-02 S at degree 2 (`20261016T070008Z_WF-10-step4-H20.log:6836`), while per phantom watt \|B₁⁺\| moved +0.40 % and SAR10g(c\*) +1.58 %. Entry 39 found the degree-2 root. The tuned figures are therefore re-read **at that root** before either reading of entry 35 stands. This window's 10 g module regenerates the 18-point lattice table that item 491's H27 parser reads, so it also serves "H27 at the degree-2 node".

**Price (measured in the family; scaling stated):**
- **Entry 35**, with the same two modules, mesh (369 387 cells), degree (2 570 024 unknowns) and `-n 16`: **1 325 s**, `memory.peak` **86.81 GiB**. Its first module took 526.8 s: build 145.2 s, node pair 367.1 s (`…H20.log:6854`).
- **Scaling:** the knob changes one capacitance in Z_RC / Z_th / V_s, and matrix size and sparsity are unchanged, so the scale is 1.0.
- **H26 on the same mesh:** degree-2 node pairs ran 380–437 s (×1.04–1.19 of 367 s), with a peak of 91.09 GiB (`20261017T070008Z_WF-10-step1e-H26.log:8641`, `:10072`).
- **Predicted 1 300–1 750 s, 87–110 GiB.** Timeout 14 400 s; 754 GiB (`xxl`) or 512 GiB (`xl`) limit.

**Readout (record).** Entry 35's gate set, asserted at degree 2:
- (A0), (iv), (Deg), (P_RC), NEG node, (P_w), (Lin), (Σ_P), (Sense), and (M_B) at `C4_COVARIANCE_BAND`;
- the 10 g module's (S_rep), against the `h9p-node-115.5209-d2.json` its sibling writes earlier in the window;
- coverage (i), the mirror (iii-a) / (iii-b) at `C4_COVARIANCE_BAND`, and containment (iv).

What it prints:
- **Primary moves:** mean \|B₁⁺\|_ccw and SAR10g(c\*), from the d1 root to the d2 root. The d1 values are 2.884091e-07 T and 4.852895959e-03 W/kg (`20261013T111006Z_WF-10-step4b-W1-regate.log:6390`). Entry 35's degree-2 values at 118.2437 pF are printed beside them: 2.210567e-07 T (`…H20.log:6850`) and 2.873209e-03 W/kg (`:13453`).
- **Per watt reaching the phantom:** the same two moves, with P_phantom/P_src beside the d1 value.
- **Y_in at the node, all four readings (printed, never asserted).** Predicted from H26's grid: |Im Y_in(P17, (1,+1j))| ≤ 3e-5 S, and the other three readings ≤ 2e-4 S in magnitude. An Im(P17, (1,+1j)) above 1e-3 S says the window is not at the root.
- **The 18-point 10 g lattice table** (for H27 at the root), and the z_c spread beside 359's and entry 35's.
- **Price:** per-module time and RSS, and `memory.peak`.

**Decision rule (pre-registered, mirroring entry 35's; the first weekly after the window reads it):**
- Every gate green and the d1-root → d2-root moves of mean \|B₁⁺\| and SAR10g(c\*) both ≤ 5 % ⇒ the human-scale tuned figures may be quoted at degree 1, with the root-to-root order move beside them.
- Either move > 5 % ⇒ they are quoted at degree 2 only, from this window, at C_in(d2).
- A gate red ⇒ a known-issues entry: the first degree-2 red at the degree-2 root on the tuned copper F-human.

No band moves.

**Command** (provisional shape; **the review that marks READY confirms the knob name and value against item 511's landing commit before copying**, and sets the profile and service to the night's tier, as entry 39's READY review did). It is entry 35's command (`:1781` of this file) with exactly three changes:
1. `XL_CHUNK` goes from `WF-10-step4-H20` to `WF-10-step4-H28`.
2. The raw-log name goes from `wf10-h20-raw.log` to `wf10-h28-raw.log`.
3. `FEM_EM_WF10_NODE=d2root` is added after `FEM_EM_WF10_STEP2_SAR=1`.

```
XL_CHUNK="WF-10-step4-H28"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf10 && R=/workspace/logs/wf10-h28-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF10_DEGREE=2 FEM_EM_WF10_STEP2=1 FEM_EM_WF10_STEP2_SAR=1 FEM_EM_WF10_NODE=d2root timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_wf10_f_human_tuned_b1_sar.py tests/validation/test_wf10_f_human_tuned_sar10g.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node record `logs/wf10/h9p-node-115.5209-d2.json` is gitignored. Entry 35's `h9p-node-118.2437-d2.json` is not overwritten, and entry 38 (`MAT-15`, 10-19) reads it.

**Flag for the weekly (not a command edit):** entry 38's module computes its own Z_RC from `C_IN_F_HUMAN_PF` (`test_mat15_f_human_tuned_heating.py:360`). P22 therefore also runs at degree 2 on the degree-1 node, and its d1 → d2 move will carry the same detuning entry 35's did.

### 42. `xl` — `MAT-15` step 6 / P22 at the degree-2 root: the tuned copper F-human's heating at degree 2, C_in(d2) = 115.5209 pF (weekly 2026-10-17, §10 chain P22; supersedes entry 38)

**Status:** **READ 2026-10-19 (weekly, §10 rulings) — first branch fires:** every `MAT-15` identity green at the root, \|Im Y_in\| ≤ 1.91e-4 S on all four readings, the root-to-root move of the ball-mean θ(t\*) per accepted watt −1.27 % (10 g θ̄ −0.20 %, SAR10g(c\*_θ) +0.26 %) ⇒ P17's thermal figures may be quoted at degree 1 with the order move beside them (`MAT-15` row); the `cem43` / six-minute clause is UNCOUNTABLE from this module (the 10-18 ruling ratified) → §10 P40 on the saved degree-2 series; the price (551 s / 86.53 GiB, OUTSIDE-below both brackets) is a record for the family's next prediction; no re-run. *Was:* **RUN `20261019T070008Z_MAT-15-step6-P22-d2root.log`** (2026-10-19 02:00; Status 0, 551 s, 86.53 GiB; ledger row filled by the 03:00 review; the decision rule is the 10-19 weekly's). *Was:* **QUEUED 2026-10-18 (10:30 review) for Monday 2026-10-19 02:00** — `docs/testing/xl-queue.d/20261019-MAT-15-step6-P22-d2root.env` (the command below, `XL_CHUNK` / `XL_COMMAND` matched line-for-line against this entry); entry 38's superseded `20261019-MAT-15-step6-P22.env` removed in the same commit; no 10-19 / 10-20 swap (P22b landed 04:52). Knob confirmed: `FEM_EM_WF10_NODE` / `d2root` in `test_wf10_f_human_tuned_b1_sar.py:81`, `:147`, read by `test_mat15_f_human_tuned_heating.py:210`, `:340`. *Was:* PENDING PREREQUISITE — §10 **P22b** (lane 1, no solve): **(i) MET 2026-10-18 (§9 item 515, the item-515 outcome commit; `20261018T094411Z_MAT-15-step6b-W1.log`, `20261018T094613Z_MAT-15-step6b-W2-regate.log`, `20261018T094815Z_MAT-15-step6b-W3-regate.log`; names as landed: `FEM_EM_WF10_NODE`, `d2root`, `C_IN_F_HUMAN_D2_PF`)** — the step-5 module (`tests/validation/test_mat15_f_human_tuned_heating.py`) takes its node through `FEM_EM_WF10_NODE` (H28's reader on `test_wf10_f_human_tuned_b1_sar.py`) instead of its own `C_IN_F_HUMAN_PF` at `:360`, build-only (Deg-2) 2 570 024 at `-n 2` with the knob, the unset path bitwise (the step-5 and 6a logs reproduced). **The implementer landing P22b marks (i) MET in the same commit and does not write the command**; the review that sees it READY replaces `xl-queue.d/20261019-MAT-15-step6-P22.env` with the command below (`cmp` against this entry). If P22b is not on `main` by the 10-18 18:00 review, that review swaps the 10-19 and 10-20 queue files (entry 37 runs Monday 10-19, this entry Tuesday 10-20) — a night is not spent on the detuned node.

**Why:** entry 38's as queued; plus entry 35's reading — at the degree-1 node the degree-2 coil is detuned (Im Y_in +1.346e-02 S) and P_phantom/P_src falls 40 %, so a thermal column there would measure the matching term, not the heating at the production convention.

**Price (entry 38's, unchanged — the knob changes one capacitance):** predicted 600–1 200 s, 90–110 GiB (entry 35 ran 1 325 s / 86.81 GiB with two modules; this is one module plus the CN solve). Timeout 14 400 s; 512 GiB limit.

**Readout (record; `MAT-15` step 5's identities asserted at their imported bands — the energy identity, the calorimetric slope, the mirror on the 10 g source):** the ball-mean θ(t\*), the 10 g ball-mean θ(t\*) at c\*_θ and SAR10g(c\*) per accepted watt at the degree-2 root, beside P17's degree-1 values at the degree-1 root (root-to-root moves) and beside entry 35's per-phantom-watt reading; `cem43` and the six-minute maximum on the degree-2 series beside item 380's; Y_in at the node printed (an Im above 1e-3 S says the window is not at the root).

**Decision rule (pre-registered; the first weekly after the window reads it):** every identity green and the root-to-root move of the ball-mean θ(t\*) per accepted watt ≤ 5 % ⇒ P17's thermal figures may be quoted at degree 1 with the order move beside them; > 5 % ⇒ degree 2 only; an identity red ⇒ a known-issues entry. No band moves.

**Command** (entry 38's with three changes: `XL_CHUNK`, the raw-log name, and `FEM_EM_WF10_NODE=d2root` after `FEM_EM_WF10_DEGREE=2`; the READY review confirms the knob name against item 511's landing and P22b's):

```
XL_CHUNK="MAT-15-step6-P22-d2root"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/wf10 /workspace/logs/mat15_s5 && R=/workspace/logs/mat15-p22-d2root-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_WF10_DEGREE=2 FEM_EM_WF10_NODE=d2root FEM_EM_MAT15_STEP5=1 FEM_EM_MAT15_STEP5_DEGREE=2 FEM_EM_MAT15_STEP5_SERIES=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_mat15_f_human_tuned_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the degree-2 series goes to the gitignored `logs/mat15_s5/`.

### 43. `xl` — `IMP-3` step 10 / I30: the rod-tip ratio with the near-T1 edge inside the 4 mm cap ball, r/4 grain, degree 2, `-n 16` (weekly 2026-10-17, §10 chain I30)

**Status:** **READ 2026-10-24 (weekly, §10 rulings) — first branch: η_rod(T1) 3.539570 inside [3.48, 3.56] with every I14 identity green ⇒ `IMP-3`'s row carries the four values (3.482328 / 3.556262 / 3.530461 / 3.539570) as its record and known-issues 2026-10-14 retires by its own clause. Caveats recorded: the interval was derived from the three rungs it brackets, so "inside" separates little (§10 epitaph); at T1 the b4 mesh is finer than r/8's; the 1 g ball mass and the degree-2 dof match were printed, not asserted, on this path (this entry's Readout overstated them). Reading of record: the pointwise tip ratio at degree 2 is 3.48–3.56 across four meshes, a record, not a converged value; no finer rung. Price INSIDE both brackets.** *Was:* **RUN `20261023T070008Z_IMP-3-step10-I30.log`** (2026-10-23 02:00; Status 0, 13 passed ×16, 607 s, 89.75 GiB; η_rod(T1) 3.539570 inside [3.48, 3.56], every identity green ⇒ the rule's first branch; the 1 g ball mass is printed, not asserted, on this path — ledger row 2026-10-23; the 10-24 weekly reads it). *Was:* **QUEUED 2026-10-18 (10:30 review) for Friday 2026-10-23 02:00** — `docs/testing/xl-queue.d/20261023-IMP-3-step10-I30.env` (the command below, matched line-for-line). Knob confirmed against item 532's landing `2f57ad2b`: `KNOBS = ("r4s2", "r8", "r8build", "b4", "b4build")` (`test_imp3_tip_ladder_1g.py:97`), `b4` dispatched at `:711` inside the `ROD_CELL_GUARD_XL` patch. Trailing-7-day count on 10-23: 10-18 … 10-23 = six. *Was:* **(i′) MET 2026-10-18 (§9 item 532, the `IMP-3` I30a commit on `main` that carries this line; `20261018T110240Z_IMP-3-I30a-W0-import.log`, `20261018T110251Z_IMP-3-I30a-W1-b4build.log`, `20261018T110425Z_IMP-3-I30a-W2-r8build-regate.log`).** The selector `FEM_EM_IMP3_STEP7=b4` has landed in `tests/validation/test_imp3_tip_ladder_1g.py` (dispatch inside the `ROD_CELL_GUARD_XL` patch; `b4build` measured 351 478 cells exact, (T1in) green). Not queued by the implementer; the next review confirms the knob and queues. *Was:* **PENDING PREREQUISITE (2026-10-18 03:00 review) — not queued: the knob does not exist.** The command sets `FEM_EM_IMP3_STEP7=b4`, but `tests/validation/test_imp3_tip_ladder_1g.py` accepts only `KNOBS = ("r4s2", "r8", "r8build")` (`:84`) and skips otherwise (`:458–459`) — the window would run Status 0 with one skip and read nothing. Item 464 (`fee90cf8`) landed the `tip_ball_radius` knob in `src/` and the build-only `tests/mesh/test_birdcage_implant_rod_ball.py`, not the `b4` dispatch; item 465's code (the `_b4()` body) was never queued. **New prerequisite (i′): §10 I30a, §9 item 532** — `b4` and `b4build` added to `KNOBS` with a `_b4()` under the `ROD_CELL_GUARD_XL` patch (`_r8`'s shape, the patch licensed by this entry for this path only, disclosed), build-only (`b4build`) at `-n 2` asserting the 351 478 cells of item 464's W2. **The implementer landing I30a marks (i′) MET in the same commit and does not write the queue file**; the next review confirms the knob and queues this entry for Friday 10-23 (the night stays free: the trailing-7-day count on 10-23 is 10-18 … 10-23 = six). If I30a has not landed by the 10-22 03:00 review, 10-23 goes unspent or to a priced variant (step 6b rule 4), and this entry waits. *Was:* READY 2026-10-17 (weekly) for **Friday 2026-10-23 02:00** (trailing-7-day count on 10-23: 10-18 … 10-23 = six). Prerequisite (i) MET — I17 (§9 item 464, `20261015T123700Z_IMP-3-I17-W2-n2.log`): the `tip_ball_radius` knob, at 4 mm the near-T1 cells inside the ball (130 cells, every vertex ≤ 3.15 mm from the cap centre), 351 478 cells — over `ROD_CELL_GUARD_W2` (320 000), which is why this is an `xl` window and not a slot. The clerk copies the command below into `xl-queue.d/20261023-IMP-3-step10-I30.env` **after confirming the knob names (`FEM_EM_IMP3_STEP7=b4`, the 4 mm selector) against item 464's landing commit and the I30 brief in `IMP-3.md`**; the guard is patched to the `xl` value on this path only, as entry 33 did, and disclosed in the log.

**Why:** entry 33's three-rung cap-ball ladder reversed at r/8 (3.482 → 3.556 → 3.530) because the near-T1 edge sat outside the 2 mm ball; I30 reads the ratio once with the edge inside the 4 mm ball at the r/4 grain. One window, no ladder (10-14 epitaph).

**Price (measured in the family; scaling stated):** entry 33 — the r/8 rung at 342 072 cells, 2 200 500 degree-2 dofs, `-n 16`: 511 s window (build 93 s, sweep 278 s, SARs 25 s), 80.55 GiB `memory.peak`. This fixture is 351 478 cells (×1.03) at the coarser r/4 grain ⇒ **predicted 500–700 s, 75–100 GiB**. Timeout 14 400 s; 512 GiB limit.

**Readout (record; I14's identities and the 1 g ball ratio asserted at their imported bands):** η_rod(T1) at degree 2 on the 4 mm ball at the r/4 grain, printed beside the three-rung record (3.482328 / 3.556262 / 3.530461); the 1 g tip Q; the near-T1 local edge ratio; the cell count and price.

**Decision rule (pre-registered; the first weekly after the window reads it):** every identity green and η_rod(T1) inside [3.48, 3.56] ⇒ the ladder's spread is read as the ball-edge effect and the row carries the four values as its record; outside ⇒ the pointwise tip ratio is recorded as unresolved at this grain and `IMP-3`'s row says so (no finer rung, no band). An identity red ⇒ a known-issues entry.

**Command** (entry 35's durable-capture shape with the I30 brief's module and knob; provisional until the READY review confirms the names):

```
XL_CHUNK="IMP-3-step10-I30"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/imp3-i30-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_IMP3_STEP7=b4 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`.

### 44. `xxl` (Saturday 2026-10-24) — `GEO-38` C3: the untuned 64×64 of F-clinical-1.5T at degree 1, 64 MHz, `-n 16` (weekly 2026-10-17, §10 chain C3)

**Status:** **READ 2026-10-24 (weekly, §10 rulings) — every assertion green (reciprocity 1.887168e-14, σ_max 0.999970196163, `TH-14`'s copper identity 2.177e-12, the 64×64 saved); the band question ruled: the readout's `ADJACENT_SPREAD_BAND` (0.5 %) is the band of record — a disclosed substitution is a question, not a ruling (§10 epitaph) — so the "other ring, one step" class at 0.9262 % (33 others ≤ 0.3217 %) fires decision-rule branch 2: the fixture's symmetry is the finding (known-issues 2026-10-24 opened), no band moves, the module's 5 % assertion stands as the module's gate only. C4 (entry 46) proceeds unchanged; C13 (§10) reads the class's anatomy from the saved matrix. Price INSIDE the brackets, above the 480–900 s prediction: the ≈ 9 s per held drive of post-processing was unpriced — entries 47–50 price it.** *Was:* **RUN 2026-10-24 `20261024T070008Z_GEO-38-C3.log`** (Status 0, 1 233 s; every module assertion green; marked by the 2026-10-24 03:00 review — reading in `docs/planning/chunks/GEO-38.md` § "Review ruling 2026-10-24 03:00"; the class-spread band named in the readout below is returned to the 10-24 weekly). *Was:* **QUEUED 2026-10-22 (03:00 review) for Saturday 2026-10-24 02:00** — `docs/testing/xxl-queue.d/20261024-GEO-38-C3.env` (`XL_CHUNK` / `XL_COMMAND` matched line-for-line against this entry with `grep -F -x`, 2 of 2). *Was:* **READY 2026-10-21 (weekly) — `xxl` Saturday 2026-10-24 02:00 at `-n 16`.** Both prerequisites are MET (below); the tier is the `xxl` night because every `xl` night before the 10-24 weekly is spent or held (10-22 entry 41, 10-23 entry 43, 10-25 / 10-26 held for 45 / 46) and the tier is a ceiling, not a floor (entry 39's precedent) — the window is not split and nothing is carried over. **The 10-22 03:00 review copies the command below into `docs/testing/xxl-queue.d/20261024-GEO-38-C3.env` after a `cmp`-grade match against this entry** (knob `FEM_EM_GEO38_C3=1` confirmed against item 519's landing `5a4c04ac`, `test_geo38_f_clinical_matrix.py:46`; the save path `logs/geo38/c3-cu-s64-d1.npz`, `:51`). The trailing-7-day `xxl` count on 10-24 is one (nothing since 10-17). *Was:* PENDING PREREQUISITE — (i) **MET 2026-10-18 (§9 item 519, `5a4c04ac`; `20261018T183335Z_GEO-38-C3a-W0.log`, `20261018T183353Z_GEO-38-C3a-W1-build.log`)** — module `tests/validation/test_geo38_f_clinical_matrix.py`, knob `FEM_EM_GEO38_C3` (`build` = build-only, `1` = the full sweep). Was: (i) §10 **C3a** (lane 1 → standard, build-only): the 64×64 sweep module on F-clinical (`PORT-13`'s route with `PORT-19` reuse, a `FEM_EM_GEO38_C3` knob), asserting the port census (64) and the degree-1 unknown count build-only at `-n 2`; (ii) **entry 40's price** (Wednesday 10-21): a measured column with every assertion green. The implementer landing C3a marks (i) MET and does not write the command; the 10-22 03:00 review, reading entry 40, writes the command in entry 40's shape with C3a's module and knob at `-n 16` on the `xxl` profile and service, marks READY, and queues it into `xxl-queue.d/20261024-GEO-38-C3.env`. If either prerequisite is missing by then, **`xxl` 10-24 is not spent** (stated here).

**Why:** chain C's step 3 — the first solve on the clinical fixture: reciprocity, passivity and the C32 × mirror classes of the 64×64, and `TH-14`'s copper identity. §10 says "heavy if ≤ 590 s, else `xl`"; entry 40's own arithmetic (t_b + t_col > 590 s) puts it beyond heavy, and 63 held drives on a 1.5 M-unknown degree-1 system fit a 4–8 h window.

**Price (re-priced 2026-10-21 from entry 40's measured column at `-n 8`: t_b 317.0 s, t_col 199.0 s, t_held 2.03 s, 40.37 GiB; the module's own projection 644 s at `-n 8`, `20261021T070008Z_GEO-38-step1b-C1b.log:12146`):** at `-n 16` the build is gmsh-serial (t_b ≈ 300–350 s), the factorising column ≤ entry 40's (t_col [100, 250] s), 63 held drives at ≤ 2.03 s each ([60, 200] s) ⇒ **predicted 480–900 s; readout brackets window [450, 1 500] s, `memory.peak` [40, 90] GiB** (40.37 GiB at 8 ranks; 16 ranks duplicate ghost layers). Outside-below on time is a record, not a red. Timeout 14 400 s; 754 GiB limit. *Was:* **Price (to be read from entry 40):** t_b + t_col + 63 · t_held at `-n 16`; from entry 40's brackets (t_b [300, 650] s, t_col [228, 1 200] s, t_held [5, 60] s) ⇒ **predicted 900–5 600 s, 30–120 GiB**; the entry is re-priced by the READY review from the measured column.

**Readout (record):** the 64×64 at degree 1: reciprocity at `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, the 34 C32 × mirror classes at `ADJACENT_SPREAD_BAND` (asserted); `TH-14`'s copper identity (asserted); P_coil/P_in, the self / adjacent / opposite class means, the factorisation and held-drive times, `memory.peak` (printed).

**Decision rule (pre-registered; the first weekly after the window reads it):** every gate green ⇒ C3 lands, the matrix is saved for C10 (`ring_mode_rows`) and C4 is commissioned at this price; a class-spread red with reciprocity green ⇒ the fixture's symmetry is the finding (known-issues; the leg-box corner or the shield can), no band moves; a reciprocity or copper-identity red ⇒ a known-issues entry on the clinical-scale Leontovich wall.

**Command** (written 2026-10-21 by the weekly: entry 40's durable-capture shape with C3a's module and knob at `-n 16` on the `xxl` profile and service; the clerk copies it verbatim):

```
XL_CHUNK="GEO-38-C3"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c3-xxl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C3=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_matrix.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS` beyond the saved `.npz` under the gitignored `logs/geo38/` (C10 reads it; `FEM_EM_GEO38_C3_DRIVES` is left at its default 64).

*2026-10-21 03:00 review:* **(ii) MET** — entry 40 RUN green (`20261021T070008Z_GEO-38-step1b-C1b.log`, every assertion green, t_b 317.0 s, t_col 199.0 s, t_held 2.03 s, 40.37 GiB). The module's own projection reads C3 = t_b + t_col + 63·t_held = 644 s at `-n 8` (`:12146`), an `xl`-sized window; this entry's tier (the `xxl` night or an `xl` one) and its command are left to the 10-21 weekly and then the 10-22 03:00 review, as written above.

*2026-10-19 weekly:* PENDING **entry 40 only** (C3a MET 10-18). If entry 40 is not green by the 10-22 03:00 review, `xxl` 10-24 is **not spent** — I34 (entry 45) cannot stand in, since it is priced by the same entry 40.

### 45. `xl` (Sunday 2026-10-25) — `GEO-36` step 8 / I34: one degree-1 copper column on the rod-in-F-clinical ASTM box, in-memory build, `-n 8` (weekly 2026-10-19, §10 chain I34; **cost probe**, entry 40's family)

**Status:** **READ 2026-10-26 (weekly, §10 rulings) — the first branch: the rod-in-clinical-box fixture is priced (t_b 458.2 s, t_col 277.8 s, t_held 2.06 s, 41.76 GiB at `-n 8`; measured exponent 1.643 against the 0.84 assumed on 10-19, a record for the family's next prediction) and I36 is commissioned on the `xl` tier — its body is §10 I36b (a slot step: the `=1` branch's node knob, the with-rod arm and the tag-403 → gel twin, (P) / NEG / (R1) / (P_RC) at imported bands, the tip ratio printed), its window entry 48 re-dated to Friday 10-30. Wording corrected: the Readout below lists "the box containment" as asserted; (G3) lives only in the skipped build test and was not asserted on this column path (the build test asserted it on 10-19 / 10-21 on the same generator) — no re-run.** *Was:* **RUN 2026-10-25 `20261025T070008Z_GEO-36-step8-I34.log`** — 9 passed and 5 skipped (by design) on all 8 ranks, Status 0, 773 s, `memory.peak` 41.76 GiB. t_b 458.2 s, t_col 277.8 s, t_held 2.06 s, exponent 1.643, every axis inside the brackets. (P) 4.0e-13 / 1.4e-13, P_surf > 0, NEG 1.6237e-04 and (R1) 1.858e-14 are asserted, along with (A0), (D-rod), census 64 and U₁. **Box containment is not asserted on this path:** (G3) lives in the skipped build test, so the Readout's list overstates it. Read literally, this is the decision rule's first branch. Marked by the 2026-10-25 03:00 review from a `log-pathologist` reading (ledger row). The 10-26 weekly reads it, and entry 48's price is filled by the 10-26 03:00 review as written there. *Was:* **QUEUED 2026-10-22 (03:00 review) for Sunday 2026-10-25 02:00** — `docs/testing/xl-queue.d/20261025-GEO-36-step8-I34.env` (`XL_CHUNK` / `XL_COMMAND` matched line-for-line against this entry with `grep -F -x`, 2 of 2; trailing-7-day count on 10-25: 10-20, 21, 22, 23, 25 = five). *Was:* **READY 2026-10-21 (weekly) — `xl` Sunday 2026-10-25 02:00 at `-n 8`.** Both prerequisites are MET (below). **The 10-22 03:00 review copies the command below into `docs/testing/xl-queue.d/20261025-GEO-36-step8-I34.env` after a `cmp`-grade match against this entry** (knob `FEM_EM_GEO36_S8_COLUMN=1` confirmed against item 575's landing `7d065e0f`, `test_geo36_rod_f_clinical.py:406`, `:419`; the build test's own knob `FEM_EM_GEO36_ROD_BOX` stays unset — the column route builds in memory). Trailing-7-day count on 10-25: 10-20, 21, 22, 23, 25 = five. **Tighter prediction from the measured column (brackets below unchanged, pre-registered 10-19):** t_b ≈ 457 s (I35's W1b, gmsh-serial), t_col ≈ 199 × (1 606 744 / 1 525 193)^1.433 ≈ 215 s, held ≈ 15 s ⇒ window ≈ 650–750 s, `memory.peak` ≈ 40–50 GiB. *Was:* PENDING PREREQUISITE — **(ii) MET 2026-10-21** (entry 40 RUN green, read by the 03:00 review: t_b 317.0 s, t_col 199.0 s at U₁ 1 525 193, 40.37 GiB; brackets below scale from it); **(i) MET 2026-10-21 (§9 item 575, `7d065e0f`; `20261021T170043Z_GEO-36-I35-W0b.log`, `20261021T170059Z_GEO-36-I35-W1b-build-n2.log`; module `tests/validation/test_geo36_rod_f_clinical.py`, knob `FEM_EM_GEO36_S8_COLUMN` = `build` / `1`)** — W1b at `-n 2`, loadavg 0.36: t_build 457.5 s (gmsh 426.92 s), (A0) rel 0, U₁ 1 606 744, census 64, summed RSS 2.089 GiB. Was: (i) §10 **I35** (lane 1, build-only): a column knob (`FEM_EM_GEO36_S8_COLUMN`) on `GEO-36` step 8's module (`tests/validation/test_geo36_rod_f_clinical.py`, item 530) in `test_geo38_f_clinical_column.py`'s shape — one degree-1 copper column at a feed sheet with the rod in, `PRICE` printed before the post-column asserts — build-only green at `-n 2` asserting I33's 1 217 025 cells; the implementer landing it marks (i) MET in the same commit and does not write the command; (ii) **entry 40's price** (Wednesday 10-21): a measured column with every assertion green. The 10-22 03:00 review, reading entry 40, writes the command in entry 40's shape with I35's module and knob at `-n 8` on the `xl` profile and service, marks READY, and queues it into `xl-queue.d/20261025-GEO-36-step8-I34.env` (trailing-7-day count on 10-25: 10-20, 21, 22, 23, 25 = five, plus 10-26 makes six on that night). If either prerequisite is missing by then, the night goes unspent or to a priced-family variant under the daily's licence (step 6b rule 4), and this entry waits.

**Why:** chain I's first solve on a clinical fixture with an implant — the F2182-style rod in the ASTM box inside F-clinical-1.5T (I33: identities green at `-n 2 / 4`, 1 217 025 cells, gmsh 418.82 s under the 450 s rule; `20261019T…_GEO-36-step8` logs on item 530). Its price sets the tier of the rod's two-port drive (I36, unenumerated until this reads) and of the disc twin (I32). Unpriced at every tier: the first window is a cost probe, in the same family as entry 40 (the same generator, the box, one more body).

**Price (predicted from entry 40's brackets; the readout is the price):** I33's mesh is 1.06× entry 40's cells (1 217 025 against 1 147 366) and the rod adds conductor surface; build t_b entry 40's × 1.2 (gmsh 419 s at `-n 4` measured), column t_col entry 40's × (U₁ ratio)², held drive ≈ 15 s. Brackets: t_b [350, 800] s; t_col [250, 1 500] s; window [700, 2 800] s; `memory.peak` [35, 140] GiB. Timeout 14 400 s; 512 GiB limit.

**Readout (record; fixed brackets, no gate, no band):** RSS per rank and summed after the column and at the end; the measured exponent against entry 40's; the PROJECTIONS line. **Asserted, by I33's module as written plus the column's own:** before the column (A0) cells 1 217 025 within `CELL_COUNT_BAND`, the rod's conductor identities (D-rod), the box containment, the port census 64; after the column (P) ≤ 1e-6 on the driven sheet, P_surf > 0, NEG, (R1) ≤ `RECIPROCITY_BAND` on the held pair. PRICE prints before the post-column asserts.

**Decision rule (pre-registered; the 10-28 weekly reads it):** column measured and every assertion green ⇒ the rod-in-clinical-box fixture is priced and I36 (the two-port drive with / without the rod, the tip 1 g ratio at clinical scale) is commissioned at this price on the tier it implies; (P), (R1) or NEG red with the identities green ⇒ the price stands as a record and a known-issues entry on the rod's conductor model at clinical scale; (A0) or an identity red ⇒ no price, the step returns to the weekly, no second window; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly.

**Command** (written 2026-10-21 by the weekly: entry 40's durable-capture shape with I35's module and knob at `-n 8`; the clerk copies it verbatim):

```
XL_CHUNK="GEO-36-step8-I34"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && R=/workspace/logs/geo36-i34-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_COLUMN=1 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the `=1` route writes no file.

### 46. `xl` (Monday 2026-10-26) — `GEO-38` C4: the degree-1 in-model Larmor root on F-clinical-1.5T with r = 0.2 Ω, `-n 16` (weekly 2026-10-19, §10 chain C4)

**Status:** **READ 2026-10-26 (weekly, §10 rulings) — the first branch: C_in(F-clinical-1.5T, d1) exists, one − → + crossing of Im Y_in(P33, (1, +1j)) between 192.5000 and 194.0000 pF, interpolated 193.517801 pF, every gate green; C5 / C6 are posed at it. The three returned points: (i) the coarsest interval — entry 51 (10-28) is kept and refines it to 0.25 pF, and entry 50's degree-2 node is entry 51's refined crossing (fallback 193.5178 pF); (ii) the root lies 0.77 % below the circuit layer's `C_mode` where F-human's lay 0.106 % above — the sign of the circuit-vs-field offset is fixture-dependent, recorded, not posed; (iii) the counter-handed series positive at every node — expected for a quadrature feed read on its counter-rotating combination, a printed record.** *Was:* **RUN 2026-10-26 `20261026T070008Z_GEO-38-C4.log`** — 5 passed and 8 skipped (by design: the numpy coding tests skip when `FEM_EM_GEO38_C4` is set, `test_geo38_f_clinical_esr_root.py:689–690`) on all 16 ranks, `[capture] rc=0` (`:14212`), Status 0, 2 215 s, `memory.peak` 49.07 GiB (`:14211`), orphans 0 / 0. Every axis is inside the readout brackets. **Asserted green:** (A0) rel 0 and npz cells equal (`:12142`, `:12149`), census 64 (`:12143`), U₁ 1 525 193 exact (`:12146`), (P50) 1.413e-12 (`:12363`), NEG 162.3× (`:12364`, on the Z50 pair only, not per node), (Z50) 2.322787e-13 (`:12367`), and (P_RC) ≤ 5.462e-13 at all eight nodes. **Im Y_in (P33, a = (1, +1j))** reads −8.326607e-03 S at 192.5000 pF and +3.944861e-03 S at 194.0000 pF, so there is **one − → + crossing, interpolated 193.517801 pF** (`:14152`). That is 0.77 % below C10's 195.0265 pF, the opposite side from F-human's +0.106 %. The mirror series (P41, (1, −1j)) crosses in the same interval at 193.647735 pF (`:14155`). The two counter-handed series are positive at every node (`:14153–14154`). Price: t_b 357.7 s, Z50 pair 203.3 s, node pairs 196–215 s (`:14156`). **Read literally, the decision rule's first branch fires.** Three points are returned to the 10-26 weekly: the crossing sits in the grid's coarsest low interval (1.5 pF, with one node below it); it lies below C10; and the counter-handed series do not cross. Marked by the 2026-10-26 03:00 review from a `log-pathologist` reading (ledger row), with `:12142–12149`, `:12363–12367` and `:14152–14156` re-read. A 0.25 pF refinement window is entry 51 (daily licence, variant). *Was:* **QUEUED 2026-10-25 (03:00 review) for Monday 2026-10-26 02:00** — `docs/testing/xl-queue.d/20261026-GEO-38-C4.env` (`XL_CHUNK` / `XL_COMMAND` matched line-for-line against the command below with `grep -F -x`, 2 of 2). Trailing-7-day count before the window, by `bash_guard.runs_in_trailing_week` (cutoff today − 6 d): 10-20, 21, 22, 23, 25 = five, six with it. **READY 2026-10-25 (03:00 review):** (ii) is MET on both halves. C3 was read by the 10-24 weekly: every gate green, and the class-spread finding does not read the two-port feed, so "C4 proceeds unchanged". C10's 195.0265 pF is the grid centre. The review wrote the grid, price and command below from C4a's module as landed. Knobs: `FEM_EM_GEO38_C4=1`. `FEM_EM_GEO38_C4_NODES_PF` is comma-separated, sorted, distinct and round-trips at 4 decimals (`tests/validation/test_geo38_f_clinical_esr_root.py:108–125`). `FEM_EM_GEO38_C4_DEGREE` is unset ⇒ degree 1 (`:98–105`). The module reads C3's saved `logs/geo38/c3-cu-s64-d1.npz` (on disk, 68 826 B) and asserts (A0) against its `cells`. *Was:* PENDING PREREQUISITE — *2026-10-24 03:00 review:* (ii)'s C3 half — entry 44 RUN, every module assertion green (`20261024T070008Z_GEO-38-C3.log`); whether its class spread reads against `ADJACENT_SPREAD_BAND` (1.85×, red) or `OPPOSITE_SPREAD_BAND` (green, as asserted) is the 10-24 weekly's. (ii)'s C10 half: **C10 MET (`20261024T093056Z_GEO-38-C10-W1.log`)** — §9 item 611, 2026-10-24: `C_mode`(k = 1) 195.0265 pF from the saved 64×64 (`logs/geo38/c10-c-mode-d1.json`), reciprocity 1.887168e-14 / σ_max 0.999970196163 re-read; not READY (the grid and command are the 10-25 03:00 review's). — (i) **MET 2026-10-20 (§9 item 561, `11424a82`; `20261020T183418Z_GEO-38-C4a-W0.log`, `20261020T183431Z_GEO-38-C4a-W1-build.log`; module `tests/validation/test_geo38_f_clinical_esr_root.py`, knobs `FEM_EM_GEO38_C4` = `build` / `1`, `FEM_EM_GEO38_C4_NODES_PF`, `FEM_EM_GEO38_C4_DEGREE`)** — W1 at `-n 2`: (A0) rel 0, census 64, U₁ 1 525 193, 2E + 2F 7 885 674, feeds P33 / P41 (Δaz 90°), 62 ESR specs; the module carries H26a's route ((Z50), (P50), NEG, (P_RC)), not this entry's readout gates (03:00 ruling: the 10-25 review reconciles). Was: (i) §10 **C4a** (lane 1, build-only): the in-model Larmor root module on F-clinical-1.5T in H26a's pattern (`tests/validation/test_wf10_f_human_copper_esr_root.py`'s route: the ESR series-RC sheet on the 64 ring gaps, the two-port feed on the bisector, a `FEM_EM_GEO38_C4` knob taking a node list, the degree knob), build-only at `-n 2` asserting the port census 64 and U₁ = 1 525 193; the implementer landing it marks (i) MET and does not write the command; (ii) **C3's gates green** (`xxl` entry 44, Saturday 10-24) and **C10** (the circuit-layer `C_mode` of the k = 1 mode from C3's saved 64×64 — the node grid's centre). The 10-25 03:00 review, reading entry 44 and C10, writes the command with a fixed 6–8-node grid (≤ 4 decimals in pF) around C10's prediction, `-n 16`, marks READY, and queues it into `xl-queue.d/20261026-GEO-38-C4.env` (trailing-7-day count on 10-26: 10-20, 21, 22, 23, 25, 26 = six). If either prerequisite is missing by then, 10-26 goes unspent or to a priced-family variant (step 6b rule 4).

**Why:** chain C's step 4 — the first tuning statement on the clinical fixture, by the route that found F-human's root (H17 at degree 1, H26 at degree 2): the in-model root is the full-wave answer and the circuit layer is a prediction to be explained (10-07 epitaph). C5 (degree 2) and C6 (the ASTM box loaded: `\|B₁⁺\|`, SAR, hotspot) are posed at this node.

**Price (to be read from entries 40 and 44; the family's scaling):** one degree-1 factorisation per node on U₁ = 1 525 193 (entry 40's t_col) plus the Z50 pair, 6–8 nodes ⇒ predicted (8–10) × t_col + t_b; memory ≈ entry 40's peak. Re-priced by the READY review from entry 40's measured column. Timeout 14 400 s; 512 GiB limit.

**Readout (reconciled 2026-10-21 by the weekly to C4a's module as landed, item 561: H26a's root-finding route asserted at its imported bands — (A0), the port census 64 and U₁, (Z50) and (P50) on the Z50 pair, NEG, (P_RC) at each node; the circuit ↔ field closure at `REDUCTION_BAND`, (P_w) and the feed-bisector mirror are C5's / C6's, posed at the found root, not this entry's):** C_in(F-clinical, d1) beside C10's circuit prediction, the sign and size of Im Y_in at each node, per-node time and RSS; the node JSON saved under the gitignored `/workspace/logs/geo38/`.

**Decision rule (pre-registered; the first weekly after the window reads it):** one − → + crossing of Im Y_in inside the grid's own resolution ⇒ C_in(F-clinical, d1) exists and C5 / C6 are posed at it; no crossing in the grid ⇒ a known-issues entry (the first root failure on the clinical fixture) and the grid is re-centred once from the measured slopes before any second window; a gate red ⇒ a known-issues entry on the clinical-scale ESR sheet. No band moves.

**Grid (written 2026-10-25 03:00 review; 8 nodes, ≤ 4 decimals in pF):** `192.5000, 194.0000, 195.0265, 195.5000, 196.0000, 197.0000, 198.5000, 201.0000`, which spans −1.30 % … +3.06 % about C10's `C_mode`(k = 1) 195.0265 pF.
- **Where the root should sit.** On F-human the in-model ESR root sat *above* the circuit `C_mode`: +0.106 % at degree 1 (118.2437 against 118.1186 pF, H17) and +0.87 % at degree 2 (H25 / H26). So the grid is dense just above the centre, at 0.47–0.5 pF resolution across 194.0–196.0 (≈ 0.25 %). It widens to 1.0–2.5 pF on the high side, which keeps a larger clinical-scale offset inside the grid, and holds one low node (192.5) so a crossing below the centre is bracketed too.
- **No measured slope.** The grid is not re-centred from one: none exists on this fixture. A no-crossing reading triggers the decision rule's one re-centre "from the measured slopes".

**Price (written 2026-10-25 03:00 review from entry 44, the same fixture at `-n 16`):** the build t_b 350.3 s. Each node pair is one factorising solve (t_col 152.3 s plus ≈ 28 s post) and one held solve (2–11 s), ≈ 185–190 s per pair. The Z50 pair plus 8 nodes is 9 pairs ⇒ **predicted ≈ 2 000–2 100 s, `memory.peak` ≈ 40–55 GiB** (entry 44: 46.85 GiB). **Readout brackets:** window [1 400, 3 600] s, `memory.peak` [30, 90] GiB. Timeout 14 400 s; 512 GiB limit.

**Command** (entry 39's durable-capture shape with C4a's module and knobs, `-n 16`, raw log `geo38-c4-raw.log`, on the `xl` profile and service as in entry 47; written 2026-10-25 by the 03:00 review; the clerk copied it verbatim):

```
XL_CHUNK="GEO-38-C4"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c4-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C4=1 FEM_EM_GEO38_C4_NODES_PF=192.5000,194.0000,195.0265,195.5000,196.0000,197.0000,198.5000,201.0000 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs (`logs/geo38/c4-node-<pF>.json`) are gitignored.

### 47. `xl` (Tuesday 2026-10-27) — `GEO-38` C14 = C11's full branch: the untuned 64×64 of F-clinical-3T at degree 1, 128 MHz, `-n 16` (weekly 2026-10-24, §10 chain C14; entry 44's family, one named knob varied — the fixture module)

**Status:** **READ 2026-10-31 (weekly; the 10-28 session was lost to the host outage) — branch 2, as annotated 10-26: every gate green at degree 1 (a record); C15 printed the "other ring, 1 step" class's absolute spread at 1.1235× the adjacent class's, so known-issues 2026-10-27 stands with its retirement condition written tonight (the 1.5 T order reading on `xxl` entry 54, read by C18; then C5e). Entry 49 runs as queued. Full text PROJECT_PLAN §10.** *Was:* **RUN `20261027T070008Z_GEO-38-C11.log` (2026-10-27 02:00; 1353 s, 36.23 GiB, Status 0; every gate green; 12 of 34 classes above 0.5 % — decision rule branch 2, read by the 2026-10-27 03:00 review from a `log-pathologist` reading; C15 is §9 item 655, the weekly reads the entry).** *Was:* **QUEUED 2026-10-25 (03:00 review) for Tuesday 2026-10-27 02:00** — `docs/testing/xl-queue.d/20261027-GEO-38-C11.env` (`XL_CHUNK` / `XL_COMMAND` matched line-for-line against this entry with `grep -F -x`, 2 of 2; knob `FEM_EM_GEO38_C11` and save path `logs/geo38/c11-cu-s64-d1-128mhz.npz` confirmed at `tests/validation/test_geo38_f_clinical_3t_matrix.py:51–53`; trailing-7-day count on 10-27 by `bash_guard.runs_in_trailing_week` (cutoff today − 6 d): 10-21, 22, 23, 25, 26 charged = five before the window, six with it — the weekly's "10-22 … 27 = five" omitted 10-21; the budget holds either way). *Was:* **READY 2026-10-24 (weekly) — `xl` Tuesday 2026-10-27 02:00 at `-n 16`.** Prerequisites MET: C11 (the 3 T 64×64 module, build-only, landed 10-21 — `tests/validation/test_geo38_f_clinical_3t_matrix.py`, knob `FEM_EM_GEO38_C11` = `build` / `1`, save path `logs/geo38/c11-cu-s64-d1-128mhz.npz`) and C12 (§9 item 595, `20261022T183301Z_GEO-38-C12-W1-build.log`: 1 057 895 cells, census 64, U₁ 1 425 765 at 128 MHz on the provisional gel). **The 10-25 03:00 review copies the command below into `docs/testing/xl-queue.d/20261027-GEO-38-C11.env` after a `cmp`-grade match against this entry** (knob and save path confirmed against the module as landed). Trailing-7-day count on 10-27: 10-22, 23, 25, 26, 27 = five. The gel (0.5 S/m / ε_r 78) is the 10-21 weekly's provisional ruling (W3) and is disclosed in the record, not asserted as the standard's value.

**Why:** chain C's first solve on the 3 T twin — reciprocity, passivity, the C32 × mirror classes and `TH-14`'s copper identity at 128 MHz; the saved 64×64 is C15's input (the 3 T `C_mode`) and entry 49's grid centre. The 1.5 T twin's 64×64 (entry 44) ran green on 10-24.

**Price (measured in the family; scaling stated):** entry 44 — the same module shape on 1 147 366 cells / U₁ 1 525 193 at `-n 16`: t_b 350.3 s, t_col 152.3 s (+ 27.6 s post), held cadence 11.04 s × 63 (≈ 9 s of it post-processing), window 1 226 s, `memory.peak` 46.85 GiB. This fixture is 0.922× the cells and 0.935× the unknowns; at the measured degree-1 exponent 1.433 the column is ≈ 138 s, the build gmsh-serial ≈ 330 s, the held cadence unchanged ⇒ **predicted 1 050–1 250 s, 40–50 GiB; readout brackets window [450, 1 500] s, `memory.peak` [30, 90] GiB** (entry 44's). Timeout 14 400 s; 512 GiB limit.

**Readout (record):** the 64×64 at degree 1, 128 MHz: reciprocity at `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, the 34 C32 × mirror classes at `OPPOSITE_SPREAD_BAND` as the module asserts, **with the worst class spread printed and read against 0.5 % (`ADJACENT_SPREAD_BAND`) by the review — the 10-24 ruling on entry 44**; `TH-14`'s copper identity (asserted); P_coil/P_in, the class means, t_b / t_col / held cadence, `memory.peak` (printed); the matrix saved.

**Decision rule (pre-registered; the first weekly after the window reads it; second branch annotated 2026-10-26 by the weekly after known-issues 2026-10-24 retired on C13's reading — the branch's action is now C13's statistic, nothing else moves):** every gate green and every class ≤ 0.5 % ⇒ C14 lands clean and C15 / entry 49 follow; every gate green with a class above 0.5 % ⇒ C15 prints that class's absolute spread against the adjacent class's (C13's statistic) on the 3 T matrix — a known-issues row only if the absolute spread exceeds the adjacent class's, otherwise a small-entry record as on the 1.5 T twin — and C15 / entry 49 still follow; a reciprocity or copper-identity red ⇒ a known-issues entry on the clinical-scale Leontovich wall at 128 MHz, entry 49 not queued. No band moves.

**Command** (entry 44's durable-capture shape with C11's module and knob at `-n 16` on the `xl` profile and service; the clerk copies it verbatim):

```
XL_CHUNK="GEO-38-C11"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c11-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C11=1 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_3t_matrix.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS` beyond the saved `.npz` under the gitignored `logs/geo38/` (C15 reads it; `FEM_EM_GEO38_C11_DRIVES` left at its default 64).

### 48. `xl` (Monday 2026-11-02; was Friday 10-30, lost to the outage; was Wednesday 10-28 — re-dated 2026-10-26 by the weekly, 10-28 is entry 51's) — `GEO-36` step 9 / I36: the rod-in-F-clinical ASTM box under the two-port drive, degree 1, `-n 8` (weekly 2026-10-24, §10 chain I36; entry 45's family)

**Status:** **READ 2026-11-02 (weekly) — the first branch, with the pathologist's correction: every gate that ran is green on both arms ((R1) 2.618e-14 / 3.346e-15, (P) ≤ 1.823e-12, NEG ≥ 19× the band, (P_RC) 1.809e-12, (A0) / (A0m) / census 64 / U₁ / (D-rod) / (M2) / (Z) exact); box containment was not asserted in this window (the generator's refusal and (G2) are the proxies; carried onto entry 55's readout as "not asserted"); (A0m) cannot fail on this route and (P_RC) repeats (P) on one column — recorded as such. I36 lands as a record: tip 1 g Q 4.460568 / 4.428040 / 4.413534 / 4.341688 on the four balls, η_rod(T1) 3.642324, the rod's share of the (1, +1j) volume power 3.740924e-05, Im Y_in −1.255e-03 / +2.957e-03 S on P33 / P41 at the node; arm M's and arm R's S agree to ≈ 1e-7. I32 is posed at this price; entry 55 kept as queued. Full text PROJECT_PLAN §10.** *Was:* **RUN 2026-11-02 — `20261102T080008Z_GEO-36-step9-I36.log`: Status 0, 1052 s, `memory.peak` 45 519 204 352 B = 42.39 GiB, `-n 8`, 5 passed / 12 skipped on all eight ranks, no FAILED or ERROR token (`log-pathologist`, 2026-11-02 03:00 review: the `=1` body ran and passed 8/8, `:12504`; every gate that ran is green). The 11-02 weekly reads it. One correction for that weekly, the same one the 10-26 weekly made on entry 45: **the Readout's "box containment" is not asserted in this window** (G3 lives in a module this window did not collect); the generator's silent `_box_refuse` and (G2) are the proxies. Also: (3.48, 3.56) is an interval on η_rod(T1), not on tip Q. η_rod(T1) reads 3.642324 here, tip Q 4.34–4.46; the definitions differ and no band applies.** *Was:* **QUEUED 2026-10-30 (18:00 review) for Monday 2026-11-02 02:00** — `docs/testing/xl-queue.d/20261102-GEO-36-step9-I36.env` (the command below, copied verbatim with the Write tool). The 10-29 03:00 review that was to write it never ran (host outage 10-28 ≈ 18:08 … 10-30 13:51, journal §1); the 10-30 night was lost, so the entry takes the first night after entry 49's. Node **193.5996 pF** (entry 51's crossing 193.599624 at four decimals); its printed forms were checked against the private-leak guard's reference set before queueing, clear. Trailing-7-day count on 11-02 (cutoff 10-27): 10-27, 10-28, 11-01 = three, four with it. **§9 item 695 (I36c) landed green, 2026-10-30 19:30 slot (`20261031T003859Z_GEO-36-I36c-W0b.log`): `_arm`'s superposition, volume powers, 1 g balls and T1 read ran on a small mesh with `_drive` stubbed. What this window runs for the first time narrows to `_drive` under `=1` and the second solver.** **§9 item 702 (I36d) DONE 2026-10-31, 04:30 slot: `_run_drive` ran end to end on stubs at `-n 1` / `-n 2` (`20261031T093652Z_GEO-36-I36d-W0b.log`), so this file stands and the first-run risk narrows to `_drive` and the second factorisation.** *Was:* **READY 2026-10-27 (§9 item 656, 12:00 slot) — prerequisite (iii) MET: I36b's body landed in `tests/validation/test_geo36_rod_f_clinical_drive.py` (`_run_drive`; node knob `FEM_EM_GEO36_S8_DRIVE_C_PF`, required under `=1`, ≤ 4 decimals; missing ⇒ `ValueError` before gmsh) and the build branch is green at `-n 2` with the body present (`20261027T170559Z_GEO-36-I36b-W1-build.log`: 5 passed / 10 skipped on both ranks, 502 s, t_build 498.3 s, (A0) 1 217 025 rel 0, U₁ 1 606 744). The command is not written by the slot: the 10-29 03:00 review writes it in entry 45's shape with `FEM_EM_GEO36_S8_DRIVE=1 FEM_EM_GEO36_S8_DRIVE_C_PF=<node>` (unset `FEM_EM_GEO36_S8_DRIVE_DEGREE` ⇒ degree 1), `-n 8`, and queues it for Friday 10-30. The `=1` body has not been executed anywhere; it writes `/workspace/logs/geo36/i36-node-<pF>.json`.** *Was:* **PENDING PREREQUISITE (iii) = §10 I36b (weekly 2026-10-26) — the body, a slot step: the `=1` branch takes its node from a required `FEM_EM_GEO36_S8_DRIVE_C_PF` knob (the review writes 193.5178 pF or entry 51's refinement), solves arm M (the rod in, the (1, +1j) quadrature drive as one factorisation with two held columns) and arm R (the same mesh with tag 403 mapped to gel, I41's route, the matched census asserted as (A0m)), asserts (A0) / census 64 / U₁ / (D-rod) / (P) / NEG / (R1) / (P_RC) at their imported bands and prints the tip 1 g Q, η_rod(T1), the rod's power share and Im Y_in at the node; build-only green at `-n 2` ⇒ READY.** The night is **Friday 10-30**: the 10-29 03:00 review (or the 10-29 18:00 review at the latest) writes the command in entry 45's shape with `FEM_EM_GEO36_S8_DRIVE=1 FEM_EM_GEO36_S8_DRIVE_C_PF=<node>`, `-n 8`, marks READY and queues it into `xl-queue.d/20261030-GEO-36-step9-I36.env` (trailing-7-day count on 10-30: 10-25, 26, 27, 28, 29 = five). Predicted from entry 45 plus one build-free factorisation: **≈ 1 050–1 200 s, 42–50 GiB** (brackets below unchanged). If I36b has not landed by 10-29 18:00, the 10-30 night goes to a priced-family variant under the daily's licence or unspent, and this entry waits. *Was:* **PENDING PREREQUISITE (iii) — written 2026-10-26 03:00 review: the module's `=1` branch is a stub.** (ii) is MET: entry 45 ran green, and its measured price is filled below. But `FEM_EM_GEO36_S8_DRIVE=1` raises `NotImplementedError` before any build (`tests/validation/test_geo36_rod_f_clinical_drive.py:19`, `:180–181`, asserted by (K2) at `:312`). That is item 593's option (a), "I36 adds its body" (`docs/planning/chunks/GEO-36.md` § item 593, for-the-review item 1). The command this entry names would therefore end in a red before the build, so **it is not written and not queued.** **(iii)** I36's drive body has to land first, insertions only and build-only green at `-n 2`: the node, the with-rod solve and the tag-403 → gel twin, the gates of the Readout below. It is a slot step whose design (option (b)'s node knob, the twin, which (P) / NEG / (R1) apply) is the 10-26 weekly's. The 10-28 night is filled by entry 51 instead. *Was:* **PENDING PREREQUISITE — entry 45's first branch (Sunday 10-25).** (i) MET 2026-10-22 (§9 item 593, I36a: `tests/validation/test_geo36_rod_f_clinical_drive.py`, knob `FEM_EM_GEO36_S8_DRIVE` = `build` / `1`, degree knob `FEM_EM_GEO36_S8_DRIVE_DEGREE` unset ⇒ degree 1; build-only green at `-n 2`, 1 217 025 cells, U₁ 1 606 744; `20261022T170259Z_GEO-36-I36a-W1-build-n2.log`). (ii) **entry 45's measured column** (t_b, t_col, t_held, `memory.peak` at `-n 8`) with every assertion green. **The 10-26 03:00 review, reading entry 45, writes the measured price into this entry (brackets below unchanged), marks READY and queues it into `xl-queue.d/20261028-GEO-36-step9-I36.env`** (trailing-7-day count on 10-28: 10-23, 25, 26, 27, 28 = five). If entry 45 is red on (P) / (R1) / NEG with the identities green, the price stands but this entry is **not queued** (the conductor model is the finding first); if entry 45 is killed or OOM, this entry waits for the 10-26 weekly's re-pose. The 10-26 weekly may re-tier it.

**Why:** chain I's first implant solve under a drive at clinical scale — the F2182-style rod in the ASTM box inside F-clinical-1.5T under the two-port feed on the bisector (P33 / P41), identities and a printed tip ratio; its price sets I32 (the disc twin) and the tier of any clinical-scale implant heating step. No published comparand until (7c).

**Price (predicted from entry 45's brackets, to be replaced by its measurement):** the drive module solves one factorisation (two feed columns on it) at U₁ 1 606 744 ⇒ t_b [350, 800] s, t_col [250, 1 500] s, two held drives ≤ 2 × 15 s; **window [700, 2 800] s, `memory.peak` [35, 140] GiB.** Timeout 14 400 s; 512 GiB limit.

**Measured price (written 2026-10-26 03:00 review from entry 45, the same mesh and `-n 8`; brackets above unchanged):** t_b 458.2 s, t_col 277.8 s (one factorising column), t_held 2.06 s, `memory.peak` 41.76 GiB (`20261025T070008Z_GEO-36-step8-I34.log:12316`, `:12382`). One factorisation with two feed columns on it is ≈ 458 + 278 + 2 × 2 s plus the drive's post-processing. With entry 46's per-node post-processing (≈ 7 s of held solve and ≈ 10–25 s beyond the factorisation, `20261026T070008Z_GEO-38-C4.log:14156`), that is **predicted ≈ 760–850 s and 42–50 GiB for the with-rod solve alone.** A with / without-rod twin, if the weekly writes one, adds a build-free second factorisation (≈ 280–300 s). These are predictions; I36's body sets the shape.

**Readout (record; the module's gates asserted at their imported bands — (A0) cells within `CELL_COUNT_BAND`, the rod's conductor identities (D-rod), the box containment, the port census 64; (P) ≤ 1e-6 on each driven sheet, P_surf > 0, NEG, (R1) ≤ `RECIPROCITY_BAND` on the held pair):** the tip 1 g ball ratio and η_rod(T1) at degree 1 printed beside `IMP-3`'s F-small record (3.48–3.56); the rod's share of the dissipated power; per-phase time and RSS.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ I36 lands as a record (no band on the ratio), I32 is posed at this price and the clinical implant-heating step is tiered from it; (P) / (R1) / NEG red with the identities green ⇒ a known-issues entry on the rod's conductor model at clinical scale, price as a record; (A0) or an identity red ⇒ the step returns to the weekly, no second window.

**Command** (entry 45's durable-capture shape with I36a's module and knobs, written 2026-10-30 18:00 review; the clerk copied it verbatim):

```
XL_CHUNK="GEO-36-step9-I36"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo36 && R=/workspace/logs/geo36-i36-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_DRIVE=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_drive.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the JSON (`logs/geo36/i36-node-<pF>.json`) is gitignored.

### 49. `xl` (Sunday 2026-11-01; was Thursday 10-29, lost to the outage) — `GEO-38` C16: the degree-1 in-model Larmor root on F-clinical-3T with r = 0.2 Ω at 128 MHz, `-n 16` (weekly 2026-10-24, §10 chain C16; entry 46's family)

**Status:** **READ 2026-11-02 (weekly) — the first branch: C_in(F-clinical-3T, d1) exists, 44.0615 pF by the bracketing pair on P33 (1, +1j) at the grid's 0.28 pF resolution (44.0505–44.0897 across the four series), every gate green; an internal tuning statement at degree 1, no absolute or external claim. The 3 T twin's C5 / C6 are posed at it (C21's body tonight; entry 58 the window, at entry 57's refined node). The sign reading — `C_mode` 43.9719 pF 0.20 % below the root here, the 1.5 T root 0.73 % below its `C_mode` — is C26's (the reader with and without r on both twins), not a window. Entry 57 kept. Full text PROJECT_PLAN §10.** *Was:* **RUN `20261101T080008Z_GEO-38-C16.log`** (marked by the 2026-11-01 03:00 review): 5 passed, 5 skipped on all 16 ranks, Status 0, 1 738 s, 39.39 GiB; every gate green; one − → + crossing on all four series inside [43.9719, 44.2500] pF, P33 (1, +1j) 44.061525 pF by the bracketing pair (ledger row). The decision rule is the 11-02 weekly's to read; its first branch reads as firing. Entry 57 refines the bracket. *Was:* **QUEUED, re-dated 2026-10-30 (18:00 review) for Sunday 2026-11-01 02:00** — `docs/testing/xl-queue.d/20261101-GEO-38-C16.env`, the same file renamed, its text unchanged: the host was down at 02:00 on 10-29 and 10-30, so the file was never consumed. Sunday is the next `xl` night (Saturday is `xxl`). Trailing-7-day count on 11-01 (cutoff 10-26): 10-26, 10-27, 10-28 = three, four with it. *Was:* **QUEUED 2026-10-28 (03:00 review) for Thursday 2026-10-29 02:00** — `docs/testing/xl-queue.d/20261029-GEO-38-C16.env` (the command below, copied verbatim with the Write tool). Both prerequisites MET: (ii) entry 47 green (10-27, every gate, `20261027T070008Z_GEO-38-C11.log`; its 64×64 at `logs/geo38/c11-cu-s64-d1-128mhz.npz`, C12's `C11_SAVE_PATH`) and C15 DONE (§9 item 655; `c15-c-mode-d1-128mhz.json`, `C_mode` 43.9719 pF). **Grid (8 nodes, written by this review):** 42.7000, 43.2000, 43.4500, 43.6000, 43.7400, 43.9719, 44.2500, 44.7000 pF — C15's centre and seven nodes from −2.9 % to +1.7 % of it, densest where entry 51 put the 1.5 T root (0.73 % below its `C_mode`), spacing 0.34–1.14 %. Every node's printed `.12e` capacitance was checked against the private-leak guard's reference set before queueing (none informative except the centre, which is clear), because entry 51's log was held on its own node line (see entry 51). **Re-priced from entries 46 / 51 (same module family, `-n 16`):** 1.5 T node pairs 196–223 s on U₁ 1 525 193 ⇒ ≈ 0.91× on U₁ 1 425 765 ≈ 180–205 s per node; t_b 376.9 s (entry 47, this mesh); Z50 pair ≈ 175–185 s ⇒ **predicted ≈ 2 050–2 300 s, 40–50 GiB; readout brackets window [1 200, 5 400] s, `memory.peak` [25, 90] GiB.** Trailing-7-day count on 10-29: 10-23, 25, 26, 27, 28 = five, six with it. *Was:* **PENDING PREREQUISITE** — (i) MET 2026-10-22 (§9 item 595, C12: `tests/validation/test_geo38_f_clinical_3t_esr_root.py`, knobs `FEM_EM_GEO38_C12` = `build` / `1`, `FEM_EM_GEO38_C12_NODES_PF`, degree unset ⇒ 1; build-only green at `-n 2`: census 64, U₁ 1 425 765); (ii) **entry 47 green** (Tuesday 10-27) and **C15** (the 3 T `C_mode`(k = 1) from entry 47's saved 64×64 through C10's reader — §10 chain C, a lane-1 slot step the 10-27 03:00 review queues; its JSON is this entry's grid centre). **The 10-28 03:00 review, reading entry 47 and C15, writes the command with a fixed 6–8-node grid (≤ 4 decimals in pF) around C15's prediction, `-n 16`, marks READY, and queues it into `xl-queue.d/20261029-GEO-38-C16.env`** (trailing-7-day count on 10-29: 10-25, 26, 27, 28, 29 = five). If either prerequisite is missing by then, 10-29 goes unspent or to a priced-family variant (step 6b rule 4).

**Why:** chain C's step 4 on the 3 T twin — the first tuning statement at 128 MHz on a clinical fixture, by the route that found F-human's root (H17 / H26) and F-clinical-1.5T's (entry 46).

**Price (to be read from entries 46 and 47):** one degree-1 factorisation per node on U₁ 1 425 765 (0.935× entry 46's unknowns ⇒ ≈ 0.91× its per-node time at exponent 1.433) plus the Z50 pair, 6–8 nodes ⇒ predicted (6–8) × t_node + t_b; memory ≈ entry 47's peak. Re-priced by the READY review from entry 46's measured per-node time. Timeout 14 400 s; 512 GiB limit.

**Readout (H26a's root-finding route asserted at its imported bands — (A0), the port census 64 and U₁, (Z50) and (P50) on the Z50 pair, NEG, (P_RC) at each node):** C_in(F-clinical-3T, d1) beside C15's circuit prediction, the sign and size of Im Y_in at each node, per-node time and RSS; the node JSON saved under the gitignored `/workspace/logs/geo38/`.

**Decision rule (pre-registered; the first weekly after the window reads it):** one − → + crossing of Im Y_in inside the grid's own resolution ⇒ C_in(F-clinical-3T, d1) exists and the 3 T twin's C5 / C6 are posed at it; no crossing ⇒ a known-issues entry and one re-centre from the measured slopes before any second window; a gate red ⇒ a known-issues entry on the clinical-scale ESR sheet at 128 MHz. No band moves.

**Command** (entry 51's durable-capture shape with C12's module and knobs, written 2026-10-28 03:00 review; the clerk copied it verbatim):

```
XL_CHUNK="GEO-38-C16"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c16-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C12=1 FEM_EM_GEO38_C12_NODES_PF=42.7000,43.2000,43.4500,43.6000,43.7400,43.9719,44.2500,44.7000 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_3t_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs (`logs/geo38/c12-node-<pF>-128mhz.json`) are gitignored. **Not exercised before:** C12's `1` branch has never run (item 595 stopped at build-only); it composes C4a's `_run_nodes`, which ran green in entries 46 and 51.

### 50. `xxl` (Saturday 2026-10-31) — `GEO-38` C5a: one degree-2 node on F-clinical-1.5T — the degree-2 cost probe, `-n 16` (weekly 2026-10-24, §10 chain C5a; **cost probe**, priced from entry 32's F-human degree-2 column and entry 40's exponent)

**Status:** **READ 2026-10-31 (weekly) — the first branch: the clinical degree-2 column is priced. Tier ruled (the rule's letter plus the memory clause it did not write): a root grid runs on `xl` at ≤ 3 nodes per window (entry 53's shape); an OOM on `xl` moves the family to `xxl` at the same ≤ 3 nodes with no `xl` re-run; the column set (entry 54) and any root-node column set (C19) are `xxl`. Entries 53 and 54 kept as written. (Z50) is not an asserted gate at degree 2 — the Readout below is read against the module's asserts. Full text PROJECT_PLAN §10.** *Was:* **RUN 2026-10-31 `20261031T070008Z_GEO-38-C5a.log`** — Status 0, 6 315 s, 472.42 GiB `memory.peak` (`xxl-ledger.md` row, 10:30 review): every asserted gate green; the clinical degree-2 column is priced. The first branch reads as firing; node pair 3 056.9 s is under 1 h, so the rule's tier reads `xl`, at 92 % of its memory limit. **The 10-31 weekly rules.** Im Y_in is +5.389409e-02 S at the degree-1 root, so on degree-1's shape the degree-2 root lies below it (entries 53, 54; this one node alone gives a sign, not a root). **(Z50) is printed only at degree 2 (module `:678`), so "(Z50)" in the Readout below is not one of the gates this window asserts, the correction entry 39 had and this entry never got; read "every gate green" against the module's asserts.** *Was:* **QUEUED 2026-10-30 (18:00 review) for Saturday 2026-10-31 02:00** — `docs/testing/xxl-queue.d/20261031-GEO-38-C5a.env`. Node **193.5996 pF** (entry 51's crossing 193.599624 at four decimals; the fallback does not apply, entry 51 ran green). The 10-29 03:00 review named for this never ran (outage, journal §1), and the night is still ahead, so this review wrote it. The node's printed forms were checked against the private-leak guard's reference set before queueing, clear. `xxl` trailing-7-day count on 10-31: nothing since 10-24. *Was:* **PENDING PREREQUISITE — re-conditioned 2026-10-26 by the weekly: the node is entry 51's refined crossing (Wednesday 10-28; the 10-28 03:00 review reads it), fallback 193.5178 pF (entry 46's interpolation) if entry 51 is red or killed. The 10-29 03:00 review writes the node into the command, marks READY and queues it into `xxl-queue.d/20261031-GEO-38-C5a.env`.** Entry 46 was read first branch (every gate green, one crossing), so the "red on a gate" clause does not fire. *Was:* **PENDING PREREQUISITE — entry 46's reading (Monday 10-26 02:00 → the 10-26 weekly / 10-27 03:00 review).** (i) MET 2026-10-20 (§9 item 561, C4a: `tests/validation/test_geo38_f_clinical_esr_root.py`, knobs `FEM_EM_GEO38_C4` = `build` / `1`, `FEM_EM_GEO38_C4_NODES_PF`, `FEM_EM_GEO38_C4_DEGREE`; build-only at `-n 2` asserted U₂ = 2E + 2F = 7 885 674). (ii) **the node:** C_in(F-clinical, d1) if entry 46 finds a crossing, else C10's 195.0265 pF — one node only. **The 10-27 03:00 review writes the node into the command, marks READY and queues it into `xxl-queue.d/20261031-GEO-38-C5a.env`** (the `xxl` trailing-7-day count on 10-31 is one: nothing since 10-24). If entry 46 is red on a gate (not merely "no crossing"), 10-31 is **not spent** on this entry.

**Why:** the operator's ladder step (5) — degree 2 on `xl` / `xxl`, priced by step 1 — has no measured degree-2 price on the clinical fixture; the 10-21 weekly's re-read (the module's 1 061 s / 231–405 GiB used F-human's exponent 0.84; at the measured 1.433 the degree-2 column is ≈ 1 774 s, `xxl`-only) is a projection. One degree-2 factorisation with the Z50 pair on it is the smallest window that measures it, and it runs on the route C5 proper will use.

**Price (measured in the family; scaling stated):** entry 32 — one degree-2 copper column on F-human (2 570 024 unknowns) at `-n 16`: 514 s window, build 122.0 s, factorising drive 287.1 s, `memory.peak` 89.31 GiB (`20261013T070008Z_WF-12-step4-H16.log`). This fixture has 7 885 674 degree-2 unknowns (3.07×); with the factorisation's time scaling N^1.5–2.0 and memory N^1.2–1.5 the drive is 1 550–2 700 s and 340–480 GiB, the build gmsh-serial ≈ 350 s (entry 44), the Z50 pair two solves on one factorisation (≤ 60 s) ⇒ **predicted 1 900–3 100 s, 340–480 GiB; readout brackets window [1 200, 7 200] s, `memory.peak` [200, 700] GiB.** Timeout 14 400 s; 754 GiB limit. **Killed or OOM ⇒ the last flushed phase line (build / assembly / factorisation) is the lower bound and the price is recorded as "above N GiB at 16 ranks"; the step returns to the weekly** (a rank-count or quarter-symmetry route is then the question, `TH-20`'s precedent).

**Readout (record; C4a's gates asserted at their imported bands at the one node — (A0), census 64, U₂ 7 885 674 exact, (Z50), (P50), NEG, (P_RC)):** Im Y_in at the node at degree 2 beside entry 46's degree-1 value at the same capacitance (the order move of the admittance, printed); build / assembly / factorisation / solve times, RSS per rank and summed, `memory.peak` — **the price is the readout.**

**Decision rule (pre-registered; the first weekly after the window reads it):** every gate green ⇒ the clinical degree-2 column is priced; C5 proper (the degree-2 root grid, then the degree-2 column set) is tiered from it — `xxl` if a node costs > 1 h, otherwise `xl`; a gate red ⇒ a known-issues entry (the first degree-2 red on the clinical fixture), price as a record; killed / OOM ⇒ as above. No band moves.

**Command** (entry 46's durable-capture shape with C4a's module and knobs on the `xxl` profile and service; the 2026-10-30 18:00 review filled the node and copied it verbatim):

```
XL_CHUNK="GEO-38-C5a"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c5a-xxl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C4=1 FEM_EM_GEO38_C4_NODES_PF=\"193.5996\" FEM_EM_GEO38_C4_DEGREE=2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSON is gitignored.

### 51. `xl` (Wednesday 2026-10-28) — `GEO-38` C4r: entry 46's root bracket refined to 0.25 pF on F-clinical-1.5T, degree 1, `-n 16` (daily 2026-10-26 03:00, step 6b rule 4 — **priced-family variant** of entry 46, one named knob varied: `FEM_EM_GEO38_C4_NODES_PF`)

**Status:** **READ 2026-10-31 (weekly) — branch 1: C_in(F-clinical-1.5T, d1) = 193.5996 pF at 0.25 pF resolution is the record; entry 46's 193.5178 sat 0.042 % low (Im Y_in convex across the old bracket), so on every root grid from tonight the bracketing pair's interpolation is the reading, never a multi-point linear fit. The 193.5996 node in entries 48 / 50 / 52 is ratified. The held log is the operator's to land. Full text PROJECT_PLAN §10.** *Was:* **RUN 2026-10-28 `20261028T070008Z_GEO-38-C4r.log` — held by the private-leak guard in the gitignored `logs/leak-held/` (a false match on one input-echo line, `:13037`; the operator rules on landing it).** Status 0, 1452 s, 48.49 GiB, every gate green; P33 (1, +1j) crosses once in [193.50, 193.75] at **193.599624 pF** ⇒ branch 1 (ledger row; `log-pathologist`, COUNTABLE). Entry 50's node is written from it on 10-29. *Was:* **KEPT 2026-10-26 (weekly): not withdrawn — the root's precision is inherited by entries 50 / C5 / C6, and the family is priced; the refined crossing is entry 50's node (entry 50's status). The decision rule's first bullet is read with that: entry 50's node is written on 10-29 from this window, not on 10-27 from entry 46.** *Was:* **QUEUED 2026-10-26 (03:00 review) for Wednesday 2026-10-28 02:00** — `docs/testing/xl-queue.d/20261028-GEO-38-C4r.env` (copied verbatim from the command below with the Write tool). **The 10-26 weekly may withdraw it** by deleting that file before 10-28 02:00; if it does, the night goes unspent. Written under the daily's variant licence because entry 48 cannot run (its `=1` branch is a stub, entry 48's status) and the 10-28 night would otherwise idle. Trailing-7-day count on 10-28 (cutoff 10-22): 10-22, 23, 25, 26, 27 = five charged, six with it.

**Why:** chain C4's question is C_in(F-clinical, d1), and entry 46 answered it on its coarsest low interval. The crossing lies between 192.5000 and 194.0000 pF (1.5 pF ≈ 0.77 %, one node below it), and the interpolated 193.517801 pF rests on two points (`20261026T070008Z_GEO-38-C4.log:14152`). Entry 50 (C5a, `xxl` 10-31) takes C_in(d1) as its degree-2 node, and C5 / C6 are posed at it. So the root's precision is what the next three windows inherit. Four nodes inside the bracket read the crossing at 0.25 pF (≈ 0.13 %) and show whether Im Y_in is linear across it. This window poses no new physics, no new fixture and no band.

**Price (measured in the family; scaling stated):** entry 46 on the same mesh, module and `-n 16`: t_b 357.7 s, Z50 pair 203.3 s, node pairs 196.3–214.5 s (`:14156`), `memory.peak` 49.07 GiB (`:14211`). Four nodes plus the Z50 pair is ≈ 358 + 203 + 4 × 205 ⇒ **predicted ≈ 1 350–1 450 s, 45–55 GiB; readout brackets window [800, 3 600] s, `memory.peak` [30, 90] GiB.** Timeout 14 400 s; 512 GiB limit.

**Grid (4 nodes, ≤ 4 decimals in pF, none equal to an entry-46 node, so no entry-46 node JSON is overwritten):** `193.0000, 193.2500, 193.5000, 193.7500`.

**Readout (record; C4a's gates asserted at their imported bands, as in entry 46: (A0), census 64, U₁ 1 525 193, (Z50), (P50), NEG on the Z50 pair, (P_RC) at each node):**
- Im Y_in at each node on all four series (P33 / P41 × (1, ±1j)).
- The module's interpolated root on the recorded series. Beside it, entry 46's 193.517801 pF, and the root from a linear fit through the six points 192.5000 … 194.0000 (entry 46's two plus these four), computed by the reading review from the node JSONs.
- Whether either counter-handed series crosses in [193.0, 193.75].
- Per-node times and `memory.peak`.

**Decision rule (pre-registered; the first weekly after the window reads it):**
- Every gate green and one crossing inside [193.0, 193.75] ⇒ C_in(F-clinical, d1) is recorded at this grid's resolution, beside entry 46's interpolation. Entry 50's node is written on 10-27 from entry 46, before this window runs. The 10-28 weekly reads the two values against each other.
- Every gate green and the crossing outside [193.0, 193.75] (that is, in [192.5, 193.0] or [193.75, 194.0]) ⇒ entry 46's interpolation is off by more than 0.25 pF. The bracket is recorded, and no second window is bought.
- A gate red where entry 46's was green ⇒ a known-issues entry on the module's reproducibility at `-n 16`.

No band moves.

**Command** (entry 46's durable-capture shape, module and knobs, with only the node list changed and the raw-log name; the clerk copies it verbatim):

```
XL_CHUNK="GEO-38-C4r"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c4r-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C4=1 FEM_EM_GEO38_C4_NODES_PF=193.0000,193.2500,193.5000,193.7500 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs (`logs/geo38/c4-node-<pF>.json`) are gitignored.

### 52. `xl` (Tuesday 2026-11-03; was Sunday 11-01, which entry 49 takes after the outage) — `GEO-38` C6: the ASTM-box-loaded F-clinical-1.5T under the two-port drive at the degree-1 root, no rod — `|B₁⁺|` and the located 10 g hotspot, `-n 8` (weekly 2026-10-26, §10 chain C6; entry 45's family — the same mesh with the rod mapped to gel)

**Status:** **RUN `20261103T080009Z_GEO-38-C6.log`** (2026-11-03 02:00, 868 s, Status 0; ledger row filled by the 2026-11-03 03:00 review; **the 11-04 weekly reads it** and poses C7 and prices entries 58 / 59). **Ratified 2026-10-31 (weekly): the sense (the target is mean |B₁⁺|_ccw, ccw = (1, −1j); the body runs the lattice under ccw and prints the (1, +1j) ball at c\* — "(1, +1j)" in the 10-26 text was the label of the series the root was read on), the 10-28 18:00 daily's 18-point lattice, and the coverage identity on tags (3, 403) with no conductor limit. The command and node stand.** **QUEUED 2026-10-31 (03:00 review) for Tuesday 2026-11-03 02:00** — `docs/testing/xl-queue.d/20261103-GEO-38-C6.env`, the command below copied verbatim (both lines matched against this file, `bash -n` clean). Trailing-7-day count on 11-03, cutoff 10-28: 10-28, 11-01, 11-02 = three, four with it. **§9 item 701 (C6a‴) DONE 2026-10-31, 04:30 slot: `_run_box` ran end to end on stubs at `-n 1` / `-n 2` (`20261031T093245Z_GEO-38-C6a-W0b-r4.log`), so this file stands and the first-run risk narrows to `_drive` and the factorisation.** *Was:* **READY (2026-10-30 19:30 slot, §9 item 694 landed: the `1` branch is composed, K1–K7 green, `20261031T003506Z_GEO-38-C6a-W0b-r3.log`).** The command below is written and final. **Two wording points for the 10-31 weekly (`GEO-38.md` "Review ruling 2026-10-30 18:00"):** the lattice runs under ccw = (1, −1j), the weekly's "(1, +1j)" printed beside it at c\*, because "(1, +1j)" and "|B₁⁺|_ccw" are opposite senses on this pair; coverage (i) runs on the gel as tags (3, 403), with no conductor limit on the ball. *Was:* **PENDING PREREQUISITE** — *(2026-10-30, §9 item 693: the module and its `build` branch landed, but its `1` branch still raises `NotImplementedError`, so (i)'s readouts are not composed and a command written today would stop at once; (ii) met. Not READY until a step composes them.)* (i) §10 **C6a** (a slot step, no solve → standard, build-only at `-n 2`): a module `tests/validation/test_geo38_f_clinical_box_loaded.py` (knob `FEM_EM_GEO38_C6` = `build` / `1`, the node knob as I36b's) that takes I36b's arm R — `GEO-36` step 8's rod-in-ASTM-box mesh (1 217 025 cells, U₁ 1 606 744) with tag 403 mapped to gel, I41's route — under the (1, +1j) quadrature drive on P33 / P41 with r = 0.2 Ω on all 64 gaps, and composes `WF-10` step 4's readouts on the box: mean \|B₁⁺\|_ccw over the box's 18-point lattice, the located 10 g hotspot (`mass_averaged_sar` on the box's balls), the coverage identity (i); (ii) **I36b landed** (the arm-R route it imports). **The first review after both land writes the command in entry 45's shape with the module and knob at `-n 8`, marks READY and queues it into `xl-queue.d/20261101-GEO-38-C6.env`** (the Sunday night; trailing-7-day count on 11-01: 10-27, 28, 29, 30 = four, five with it). If either prerequisite is missing by the 10-31 18:00 review, the 11-01 night goes to a priced-family variant under the daily's licence or unspent, and this entry waits.

**Why:** chain C's step 6 — the ASTM box loaded, `|B₁⁺|`, 1 g / 10 g SAR, a located hotspot — on the fixture the 2026-10-14 directive names for all new human-scale work, at the root entry 46 found; the first clinical-scale safety quantity, and the fixture C7 (the `ANS` case) is posed on. The rod-in-box mesh is reused because its box *is* the ASTM load (I33 placed it as the standard places it); mapping the rod to gel is I41's gated route, so no new mesh is bought.

**Price (measured in the family; scaling stated):** entry 45 on the same mesh at `-n 8`: t_b 458.2 s, one factorising column 277.8 s, held 2.06 s, 41.76 GiB; entry 41's SAR / `|B₁⁺|` post-processing on F-human (369 387 cells) 82.4 s. One factorisation with two held columns plus the box's lattice and balls ⇒ **predicted ≈ 850–1 000 s, 42–50 GiB; readout brackets window [600, 2 400] s, `memory.peak` [30, 120] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands — (A0) cells 1 217 025 within `CELL_COUNT_BAND`, (A0m) the matched census exact, census 64, U₁ exact; (P) and NEG (`TH-14`'s copper identity); (R1) on the 2×2 at `RECIPROCITY_BAND`; (P_RC) at the node at `DISCRETE_IDENTITY_RTOL`; coverage (i) power / mass at 1e-10):** mean \|B₁⁺\|_ccw over the lattice and its z_c spread; SAR10g(c\*) and its location (r, dz, az); P_box / P_src and the per-box-watt figures beside the per-accepted-watt ones (the normaliser stated, ruling 2026-10-12); Im Y_in at the node (the coil was tuned on C4's stand-in cylinder — the box's detuning is printed, never a red); per-phase time and RSS. **No mirror gate** (the box has no feed-bisector plane) and **no class gate** (the box breaks C32): the maps are records.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ C6 lands as a record — the first clinical-scale `|B₁⁺|` / 10 g figure at degree 1 at the root, an **internal statement** (no tuning, compliance, homogeneity or external claim) — and C7 (the `ANS` case on F-clinical) and P33 (its thermal half) are posed by that weekly; (P) / (R1) / NEG red with the identities green ⇒ a known-issues entry on the clinical-scale conductor model with the box load, the maps as records; (A0) / (A0m) or an identity red ⇒ the step returns to the weekly, no second window; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command** (entry 45's durable-capture shape with C6a's module and knobs, written 2026-10-30 18:00 review, node 193.5996 pF; the clerk copies it verbatim once READY):

```
XL_CHUNK="GEO-38-C6"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c6-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C6=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_box_loaded.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the maps and `logs/geo38/c6-node-<pF>.json` are gitignored.

### 53. `xl` (Wednesday 2026-11-04) — `GEO-38` C5b: three degree-2 nodes around the extrapolated root on F-clinical-1.5T, `-n 16` (daily 2026-10-31 10:30, step 6b rule 4 — **priced-family variant** of entry 46, polynomial degree and `FEM_EM_GEO38_C4_NODES_PF` varied; the measured point at the target order is entry 50, same mesh, same module, same width)

**Status:** **RUN 2026-11-04 `20261104T080008Z_GEO-38-C5b.log`** — Status 0, 13 351 s, `memory.peak` 480.47 GiB, every gate green; P33 (1, +1j) crosses once in [185.75, 187.0] at 185.763346 pF (`:13037`, a record); all four series cross (`:13038–13040`), against the rule's degree-1 expectation; first branch, the 11-04 weekly reads it (ledger row, `log-pathologist` reading, 10:30 review). **KEPT 2026-10-31 (weekly) as written — grid and rule unchanged; this is the `xl` ≤ 3-node shape the tier ruling on entry 50 names, and its successor is entry 56 (the 0.25 pF refinement, or the re-centred grid if no crossing).** **QUEUED 2026-10-31 (10:30 review) for Wednesday 2026-11-04 02:00** — `docs/testing/xl-queue.d/20261104-GEO-38-C5b.env`, the command below copied verbatim. **The 10-31 weekly may withdraw or re-grid it** by editing or deleting that file before 11-04 02:00; if it does, the night goes to the weekly's choice. Written under the daily's variant licence because `xl` was one entry short of the floor and entry 50's price is now measured. Trailing-7-day count on 11-04 (cutoff 10-29): 11-01, 11-02, 11-03 = three, four with it.

**Why:** entry 50's decision rule says C5 proper starts with the degree-2 root grid, and C6's and I36's readings at degree 2 need that root. Im Y_in at 193.5996 pF is positive on all four readings, so on degree-1's shape the root lies below the node, as on F-human (entry 39). Entry 50 gave a sign and no degree-2 slope.

**Grid (3 nodes, ≤ 4 decimals in pF): `184.5000, 185.7500, 187.0000`.** Written from `logs/geo38/c4-node-193.5996-d2.json` against the degree-1 node JSONs. Each degree-2 reading sits just above the degree-1 reading at 201.0000 pF. Linear extrapolation from the two top degree-1 nodes puts them at an equivalent degree-1 capacitance of 201.28–201.55 pF. That is a move of −7.68 to −7.95 pF (≈ −4.0 %). The translated root is 185.65–185.92 pF, and the proportional reading gives 185.96–186.21 pF. Entry 39's same one-point reading bracketed F-human's root (115.49 translated, 115.55 proportional, 115.52 found). The `log-pathologist` ruling on entry 50 (journal §4) narrows and shifts that: the degree-1 secants fall from 9.5e-3 to 5.2e-3 S/pF going up in C, so the linear reading is biased low, and a quadratic through 197, 198.5 and 201 pF gives translated roots of 185.55–185.83 pF. Re Y on the crossing series maps to ≈ 202.0 pF, which implies ≈ 185.2 pF, with loss differences unseparated. **The estimates therefore span ≈ 185.2–186.2 pF.** Here the move is ≈ 3× F-human's in pF and n = 1 is the only precedent, so this is an estimate, not a prediction with a band. The risk is lopsided: the root falls below 185 pF on a slope only 7–9 % low, and the high side needs ≈ +20 %. A first grid of 185 / 186 / 187 pF was dropped for that reason. This grid leaves ≈ 0.7 pF below and ≈ 0.8 pF above the span. A 1.25 pF spacing costs little, because the secant is good to ≈ 0.08 pF at 1.5 pF (entries 46 and 51) and a refinement follows anyway. The nodes' printed `.4f` and `.12e` forms were checked against the private-leak guard's reference set before queueing: clear.

**Why three nodes and not C5's full grid.** Entry 50's node pair is 3 056.9 s, so 3 nodes fit an `xl` window and 8 fit an `xxl` one only with ≈ 7.7 h of its 8 h. Three bracket the estimate and measure the degree-2 slope, which entry 50's one node cannot, and a refinement is then a priced variant.

**Price (measured in the family; scaling stated):** entry 50 on the same mesh, module, degree and `-n 16`: build 341.8 s, Z50 pair 2 892.1 s, node pair 3 056.9 s, `memory.peak` 472.42 GiB (`20261031T070008Z_GEO-38-C5a.log:12595`, `:12650`). One Z50 pair and three node pairs: ≈ 342 + 2 892 + 3 × 3 057 ⇒ **predicted ≈ 12 400 s (3.45 h), 472–490 GiB; readout brackets window [9 600, 14 400] s, `memory.peak` [440, 512] GiB.** Entries 46 and 51 show no memory growth with node count at degree 1 (nine factorisations 49.07 GiB, five 48.49 GiB), and the ruling confirms by three routes that the Z50 factor is freed before the node factorises. **This window has four factorisations to entry 50's two, and the budget tracks MUMPS INFOG(19), which varied with the node** (482 667 and 498 030 MB in entry 50, `:12482` for the second; 7.5 % max / min in entry 39). The 512 GiB limit leaves ≈ 540.5 GB of INFOG(19), +8.5 % over the top draw, or +5.4 % if the Z50 pair carried the peak, against a family max / mean of +3.7–3.8 %. It probably fits; that is not proven, and an OOM is covered by the decision rule. Timeout 14 400 s, the `xl` ceiling, which is a 14 % margin; factorisations spread 7 % in entry 50 (2 524.98 / 2 700.21 s, `:12285`, `:12506`). 512 GiB limit, 92 % used in entry 50.

**Runs past the 04:30 slot.** This will be the first window to overlap a slot. `OPS-54`'s row says it queues nothing until a measured window has overlapped one; a window's Python imports resolve at start, so the exposure is late imports only. The reading review says whether the overlap showed anything.

**Readout (record; C4a's gates asserted at their imported bands, as in entry 50 — (A0), census 64, U₂ 7 885 674 exact, (Deg), kinds and (P_RC) at each node, (P50) and NEG on the Z50 pair; (Z50) printed only):** Im Y_in at each node on all four series; the module's sign-change root on each series, **which the module never asserts, so a green window does not mean a root was found**; the degree-2 slope against the degree-1 slope at the root (9.78e-03 S/pF on P33 (1, +1j) over [193.50, 193.75]); per-node times and `memory.peak`.

**Decision rule (pre-registered; the first weekly after the window reads it):**
- Every gate green and one − → + crossing on P33 (1, +1j) inside [184.5, 187.0] ⇒ C_in(F-clinical, d2) is recorded at this grid's resolution, as a record, and a 0.25 pF refinement (entry 51's pattern) or the column set is posed at it. P41 (1, +1j) and P33 (1, −1j) never cross at degree 1 on any of the 12 nodes, so "none" on those two matches degree 1; F-human's "all four cross" is not the expected outcome here.
- Every gate green and no crossing on P33 (1, +1j) inside the grid ⇒ the extrapolation is off by more than ≈ 1 pF. The measured slope re-centres one more window; not a known-issues entry until that second grid misses.
- A gate red ⇒ a known-issues entry on the degree-2 fixture; the Z50 pair's green in entry 50 is the comparison.
- Killed at the 14 400 s timeout or OOM ⇒ the node JSONs already written stand as records (the module writes one per node), the price is recorded, and the step returns to the weekly.

No band moves.

**Command** (entry 50's durable-capture shape, module and knobs on the `xl` profile and service, with only the node list and the raw-log name changed; the clerk copied it verbatim):

```
XL_CHUNK="GEO-38-C5b"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c5b-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C4=1 FEM_EM_GEO38_C4_NODES_PF=184.5000,185.7500,187.0000 FEM_EM_GEO38_C4_DEGREE=2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs (`logs/geo38/c4-node-<pF>-d2.json`) are gitignored.

### 54. `xxl` (Saturday 2026-11-07) — `GEO-38` C5d: the 64×64 of F-clinical-1.5T at degree 2, 64 MHz, `-n 16` (daily 2026-10-31 10:30, step 6b rule 4 — **priced-family variant** of entry 44, polynomial degree varied; the measured point at the target order is entry 50 on the same mesh)

**Status:** **READ 2026-11-07 (weekly) — the second branch, as the 10:30 review read it: every gate green at degree 2 ((Deg) 7 885 674 exact, reciprocity 1.488791e-13, σ_max 0.999932072111, (P) ≤ 2.734e-12, (iv′) 6.4235e-03 against 0.05), one class of 34 above `ADJACENT_SPREAD_BAND` — other ring, 1 step, 0.6424 % against entry 44's 0.9262 % at degree 1. The degree-2 1.5 T 64×64 is a record (`c3-cu-s64-d2.npz`, C18's and U129's input; `C_mode`(d2) 185.4048 pF); the fixture's symmetry stays the finding; no band moves. Known-issues 2026-10-27's ≥ 2× clause does not fire (the class falls ×1.44, the median class ×1.36; C18) — the asymmetry is not degree-1 discretisation, so the 3 T degree-2 window is not commissioned on it and the next reading is §10 C35 (the r2 h-rung at degree 1). The log stays leak-held (one MUMPS timing line; dashboard item 3).** *Was:* **RUN 2026-11-07 `20261107T080008Z_GEO-38-C5d.log`** (Status 0, 5364 s; the log is held unchanged in the gitignored `logs/leak-held/` by the private-leak guard, a false match on one MUMPS timing line; the ledger row carries the reading; the 11-07 weekly reads it). *Was:* **KEPT 2026-10-31 (weekly) as written — the column set is `xxl` by the tier ruling on entry 50; its reading is also known-issues 2026-10-27's first reading (the 1.5 T classes' order move, C18) and `C_mode`(d2) beside the degree-2 root.** **QUEUED 2026-10-31 (18:00 review) for 2026-11-07, `xxl-queue.d/20261107-GEO-38-C5d.env`, command verbatim.** *Was:* READY 2026-10-31 (§9 item 705, C5c, 12:00 slot): `FEM_EM_GEO38_C3_DEGREE=2` landed on `tests/validation/test_geo38_f_clinical_matrix.py`, build-only green at `-n 2` with (Deg) 7 885 674 exact (`20261031T170146Z_GEO-38-C5c-W1-build.log:11854`). `_sweep` has never run at degree 2 through this module, so this window is its first. *Before that:* PENDING PREREQUISITE — §9 item 705 (C5c). **Read against the code by the 18:00 review (no run):** `_sweep`'s drive calls are the ones entry 50 ran at degree 2, and no defect was found. The one difference is the power partition on every column, which C4 runs only on the factorising one, so entry 50's 22.41 s excludes it. Entry 36's degree-2 post, partition included, is 4.2–5.4 s on 369 387 cells and 32 sheets; scaled ×3.1–6.2 that is ≈ 35–55 s per held column, inside this entry's 22–60 s, and the prediction stands (`GEO-38.md` "Review ruling 2026-10-31 18:00"). **The 10-31 weekly may withdraw or replace it:** the 10-26 weekly declined a degree-2 64×64 *before C5a prices degree 2*, entry 50 has now priced it, and the choice of what `xxl` 11-07 buys is the weekly's. A root-node column set on the same night needs entry 53's reading first.

**Why:** entry 50's decision rule names "the degree-2 column set" as C5's second half. One factorisation with 63 held drives gives the clinical 64×64 at degree 2 for ≈ 1.4–2.2 of entry 50's node pairs, and all capacitances are then a circuit-layer reading at no further solve. It answers three things the degree-1 matrix cannot: reciprocity and passivity at the production order on this fixture; `C_mode` at degree 2 beside entry 53's in-model root, the circuit layer being a prediction and the root the claim (10-07 epitaph); and whether the class spreads of known-issues 2026-10-24 and 2026-10-27 move with the order, which is the mesh-or-physics question the 3 T entry leaves open.

**Price (measured in the family; scaling stated):** entry 50: build 341.8 s, factorising solve 2 869.70 s, held solve 22.41 s (`:12362`). Entry 44's held cadence was 11.04 s, of which ≈ 9 s was post-processing, at degree 1. Sixty-three held drives at 22–60 s each ⇒ **predicted ≈ 4 600–7 000 s, 472–520 GiB; readout brackets window [3 600, 10 800] s, `memory.peak` [440, 700] GiB.** Timeout 14 400 s; 754 GiB limit. Entry 44 added no visible memory for the held drives (46.85 GiB against entry 46's 49.07).

**Readout (record; asserted at the module's imported bands as in entry 44 — (A0), census 64, (U1), (Deg) U₂ 7 885 674 exact, kinds 1 solve / 63 held, (P) and NEG on every column, (ii) `RECIPROCITY_BAND`, (iii) `MATRIX_PASSIVITY_CEILING`, (iv′) `OPPOSITE_SPREAD_BAND`):** the 34 class means and spreads printed beside entry 44's, with `ADJACENT_SPREAD_BAND` as the band of record for a class red (10-24 ruling); the saved `logs/geo38/c3-cu-s64-d2.npz`; held cadence and `memory.peak`.

**Decision rule (pre-registered; the first weekly after the window reads it):** every gate green and no class above `ADJACENT_SPREAD_BAND` ⇒ the degree-2 matrix lands as a record and known-issues 2026-10-27 gains the order reading; every gate green with a class above the band ⇒ the fixture's symmetry stays the finding, the class is compared with entry 44's, no band moves; a reciprocity or copper-identity red ⇒ a known-issues entry on the clinical-scale Leontovich wall at degree 2; killed or OOM ⇒ the last flushed phase line is the lower bound, and the step returns to the weekly.

**Command** (entry 44's durable-capture shape with one knob added; the clerk copies it verbatim once READY):

```
XL_CHUNK="GEO-38-C5d"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c5d-xxl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C3=1 FEM_EM_GEO38_C3_DEGREE=2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_matrix.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the `.npz` is gitignored.

### 55. `xl` (Friday 2026-11-06) — `GEO-36` step 10 / I32: the clinical implant-heating figure — arm M's SAR heats the ASTM gel box around the rod inside F-clinical-1.5T, degree 1, `-n 8` (weekly 2026-10-31, §10 chain I32; entry 48's family — the same mesh, the same drive, `MAT-15`'s thermal route added)

**Status:** **READ 2026-11-07 (weekly) — the pre-registered "slope identity red with the EM gates green" branch, as filed: the EM record stands (every EM gate, (B), (E), (E′), the global (S2) green), the two Readout corrections are ratified (the tip θ is a 1 g ball mean; Im Y_in is entry 48's print), and `rod (S2-tip)` 6.420e-03 against `S2_BAND` 1e-04 is ruled the comparand's form (`OPS-83` on the stub, `OPS-84` on I41's rod pair: the DG0 ball mean against the scheme's consistent-mass projection) — re-posed as option (b), the projection ⟨Πq⟩_B at the unmoved `S2_BAND` on one BE step at Δt₁/10, for I32 (I32b) and U127 (U127′) together. Records: tip 1 g ball-mean θ(t\*) rod 3.3567e-02 / 3.3684e-02 K, gel twin 3.0825e-02 / 3.0907e-02 K, ratio 1.0890 / 1.0898 (I41 1.10); rod power share 3.740924e-05. **Re-run as entry 63 (I32′, Friday 11-13) at this window's measured price once I32b lands.** No band moves.** *Was:* **RUN 2026-11-06 — `20261106T080008Z_GEO-36-step10-I32.log`: Status 1, 1015 s, `memory.peak` 47 368 220 672 B = 44.12 GiB, `-n 8`, orphans 0. One red: `rod (S2-tip)` 6.420e-03 against `S2_BAND` 1e-04 (`:12523`). Every EM gate, (B), (E), (E′) and the global (S2) are green; the gel twin's ungated tip reads 3.1e-03 / 1.3e-03 (`log-pathologist`, 2026-11-06 03:00 review; ledger row). This is the Decision rule's "slope identity red with the EM gates green" branch: known-issues 2026-11-06 is filed and the EM record stands. Two Readout clauses are corrected for the weekly. The tip θ is printed as a 1 g **ball mean**, not a centre value. Im Y_in is not printed by this module (entry 48 carries it). `OPS-83` (§9 item 859) measures the comparand's form on the stub. The 11-07 weekly reads this with U127.** *Was:* **KEPT 2026-11-02 (weekly): entry 48 read as its first branch (every gate that ran green; I36 a record), so prerequisite (i) is MET and this entry runs Friday as queued; the measured price below is ratified. One clause carried from entry 48: box containment is not asserted on this route (the generator's refusal and (G2) are the proxies) — read the Readout's EM gates against the module's asserts.** **QUEUED 2026-11-02 (03:00 review) for Friday 2026-11-06 02:00** — `docs/testing/xl-queue.d/20261106-GEO-36-step10-I32.env` (the command below, written with the Write tool), under this entry's own "first review after both land" clause: I32a landed 11-01, and entry 48 is green on every gate that ran (`log-pathologist`, entry 48's status). **The reading of entry 48 is the 11-02 weekly's: it may withdraw this by deleting the file, and one thing is open for it, the unasserted containment clause.** Node 193.5996 pF's printed forms were cleared against the leak guard for entry 48. Trailing-7-day count on 11-06 (cutoff 10-31): 11-01, 02, 03, 04, 05 = five, six with it, provided entry 56 runs on 11-05. *Was:* **PENDING PREREQUISITES** — (i) **entry 48 green** (Monday 11-02; its first branch: every assertion green ⇒ I36 lands as a record and I32 is posed at its price); (ii) **§10 I32a** (a slot step, no solve, smoke `-n 2`: the heating body composed in `tests/validation/test_geo36_rod_f_clinical_heating.py`, knob `FEM_EM_GEO36_S8_HEAT=1`, the node knob as I36b's, and stub-run on a small `create_box` with `_drive` stubbed — the source equals σ|E|²/ρ pointwise, the energy identity and the slope identity on a uniform source, a perturbed source named) — **(ii) MET 2026-11-01 (§9 item 707, 04:30 slot; `20261101T103929Z_GEO-36-I32a-W0b.log`); (i) still open.** **The 11-03 03:00 review (or the first review after both land, by 11-05 18:00 at the latest) writes the command in entry 48's shape with `FEM_EM_GEO36_S8_HEAT=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996`, `-n 8`, prices it from entry 48's measured window plus the thermal phase of `xl` entry 42's ledger row scaled by cells (the box 1 217 025 against F-human's 369 387), marks READY and queues it into `xl-queue.d/20261106-GEO-36-step10-I32.env`.** Trailing-7-day count on 11-06 (cutoff 10-31): 11-01, 02, 03, 04, 05 = five, six with it — the budget's last unit of the week. If entry 48 is red or killed, or I32a has not landed by 11-05 18:00, the 11-06 night goes to a priced-family variant under the daily's licence or unspent, and this entry waits for the 11-04 weekly.

**Why:** Phase 7 × Phase 8 at clinical scale — the first implant-heating figure on the production fixture, by the route `MAT-15` gated on F-small and I41 ran on the matched F-small pair (tip 1 g θ ratio 1.10): the rod's own SAR as the source on the gel it sits in, the box Robin-cooled. The closed-form thermal half is `IMP-7` (I45–I48b); the published comparand waits on W3 (7c). A record: no tissue, no compliance claim, no external check.

**Price (to be filled from entry 48; the EM half measured, the thermal half scaled):** the EM half is entry 48's window (predicted ≈ 760–850 s with-rod + ≈ 280–300 s for arm R, 42–50 GiB at `-n 8`); the thermal half is one CN transient to t\* on the gel (P1, DG0 source) — entry 42's thermal phase on F-human's phantom is the measured point, scaled linearly by cells. **Readout brackets, provisional until priced:** window [1 200, 4 000] s, `memory.peak` [35, 140] GiB. Timeout 14 400 s; 512 GiB limit.

**Measured price (written 2026-11-02 03:00 review from entry 48, the same mesh, drive and `-n 8`):** build 461.7 s, arm M 263.5 s, arm R 262.3 s, ≈ 59 s unattributed between the phases and the body's 1046.5 s, harness 1052 s, `memory.peak` 42.39 GiB (`20261102T080008Z_GEO-36-step9-I36.log:12503`, `:12626–12627`, `:12632`). **The thermal half is seconds.** Entry 42's CN read 0.55–0.84 s at `-n 16` on F-human's phantom (`20261019T070008Z_MAT-15-step6-P22-d2root.log.gz:6810–6813`), ≈ 2.8 s scaled by the entry's whole-mesh cell ratio 3.29. The closer point is I41's, 2.5 / 1.7 s on 66 693 / 51 693 sub cells (`20261025T201613Z_IMP-3-I41-W1b.log:3681`, `:3726`), against this window's ≈ 85 122 + 14 834 gel and rod cells. Two arms, each CN and BE plus the source, ball and probe reads, ≈ 30–120 s. **Predicted window ≈ 1 080–1 200 s, `memory.peak` 42–50 GiB; priced brackets [1 000, 1 800] s, [40, 70] GiB.** Timeout and limit unchanged.

**Readout (record; asserted at their imported bands — I36b's EM gates ((A0), (A0m), census 64, U₁, (D-rod), (P) / NEG on both arms, (R1), (P_RC)); the scheme's discrete energy identity ≤ `E_BAND`; the calorimetric initial slope dT/dt|₀ = SAR/c at `MAT-15`'s band on the tip 1 g ball and on the gel; (S0) census):** θ at the tip 1 g ball centre and at P37's six box probe points at t\*; the tip θ ratio arm M / arm R beside I41's F-small 1.10; the rod's power share and Im Y_in at the node (printed, I36's); per-phase time and RSS.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ I32 lands as a record (no band on the θ ratio) and the implant-heating `ANS` case (W3's comparand) is posed by that weekly if the F2182 reference is in the repo; an energy or slope identity red with the EM gates green ⇒ a known-issues entry on the thermal route at clinical scale, the EM record stands; an EM gate red ⇒ as entry 48's rule; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command** (entry 48's durable-capture shape with I32a's module and knobs, written 2026-11-02 03:00 review; the clerk copied it verbatim):

```
XL_CHUNK="GEO-36-step10-I32"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo36 && R=/workspace/logs/geo36-i32-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_HEAT=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the JSON and the fields are gitignored.

### 56. `xl` (Thursday 2026-11-05) — `GEO-38` C5b′: the degree-2 root of F-clinical-1.5T at 0.25 pF, `-n 16` (weekly 2026-10-31, §10 chain C5b′; entry 53's family — the same module, degree and width, only the node list varied)

**Status:** **READ 2026-11-07 (weekly) — the first branch: every gate green ((A0), (Deg), (P50) 1.945e-12, NEG 153.4×, (P_RC) ≤ 2.4e-12) and one − → + crossing of Im Y_in on P33 (1, +1j) — −1.555e-03 / +1.264e-03 S at 185.6250 / 185.8750 pF, 185.762904 pF by the bracketing pair (`:12814`), 0.0004 pF from entry 53's chord, −4.05 % from the degree-1 root. The rule's series is P33 (1, +1j) (entry 51's convention, every root in the family); P41 (1, −1j)'s crossing at 185.820370 pF is a record — the two-series spread 0.057 pF against 0.117 pF at degree 1 (C18). ⇒ C_in(F-clinical-1.5T, d2) = 185.7629 pF at 0.25 pF resolution (four decimals, as 193.5996), an internal tuning statement at the production order; its printed forms cleared by the leak guard (`20261108T030938Z_LEAK-CHECK-C19.log`). Entry 62 READY at it; `ANS-11` SPEC v2 carries it. Price inside both brackets (9 637 s, 482.41 GiB). No band moves.** *Was:* **RUN `20261105T080008Z_GEO-38-C5b-prime.log`** (2026-11-05 02:00 window, Status 0, 9637 s; ledger row filled by the 2026-11-05 10:30 review; **the 11-07 weekly reads it**). **RATIFIED 2026-11-04 (weekly): the 10:30 review's re-wording and its two-node straddling grid stand as queued** — the pathologist's trap is real (a 0.25 pF grid "inside the bracketing pair" would print no crossing and falsely take the second branch), the target and the decision rule are unchanged, and the first branch's root is `ANS-11`'s node (SPEC v2, the 11-07 weekly) and entry 62's. **QUEUED 2026-11-04 (10:30 review) for Thursday 2026-11-05 02:00** — `docs/testing/xl-queue.d/20261105-GEO-38-C5b-prime.env`, the command below written with the Write tool; against entry 53's consumed file only `XL_CHUNK`, the raw-log name and the node list differ (diff checked). **Grid (2 nodes): `185.6250, 185.8750`**, a 0.25 pF bracket centred on entry 53's crossing 185.763346 pF (0.138 pF above the lower node, 0.112 below the upper), none equal to an entry-53 node. Node forms cleared by the leak guard: `20261104T163729Z_LEAK-CHECK-C5BP.log`, Status 0 (an earlier grid, 186.0000 / 186.2500, cleared at `20261104T163646Z_LEAK-CHECK-C5BP.log` and was dropped, below). **Re-worded by the 10:30 review, for the weekly to ratify or revert (it may re-grid or delete the file before 02:00):** the clause "inside the bracketing pair" cannot be met with a crossing in this window's own grid. The root sits 0.0133 pF above entry 53's 185.75 node, so every 0.25 pF node inside [185.75, 187.0] that is not 185.75 reads positive. The module finds sign changes only among its own nodes, so such a grid would print no crossing and falsely take this entry's second branch (the `log-pathologist`'s trap, 10:30 journal §6b). The grid therefore straddles the crossing, inside entry 53's grid [184.5, 187.0]. The target (the root at 0.25 pF by the bracketing pair) and the decision rule are unchanged. **Two nodes, not three, on price:** entry 53's node pairs grew 2 970.4 → 3 258.1 → 3 582.0 s (`:13041`), in the MUMPS factorisation alone (2 595.6 → 2 877.6 → 3 199.8 s, `:12488`, `:12711`, `:12934`) with operations +2.1 %. Three nodes at the top figure sum to ≈ 14 264 s against the 14 400 s timeout; two give ≈ 9 460–10 690 s. *Was:* **PENDING PREREQUISITE — entry 53's reading (Wednesday 11-04 02:00; the 11-04 03:00 review reads it).** If 53 finds one − → + crossing of Im Y_in on P33 (1, +1j) inside [184.5, 187.0] pF: **the 11-04 03:00 review writes ≤ 3 nodes at 0.25 pF spacing inside the bracketing pair** (none equal to an entry-53 node), checks each node's printed forms against the leak guard, marks READY and queues it into `xl-queue.d/20261105-GEO-38-C5b-prime.env`. If 53 finds no crossing: the same review writes 53's **re-centred grid** from the measured degree-2 slope (53's own second branch, ≤ 3 nodes) into this entry instead, and the 0.25 pF refinement waits for the window after. If 53 is red on a gate, or killed / OOM: this entry is **not queued** (the tier ruling on entry 50: an OOM moves the family to `xxl` at ≤ 3 nodes, the weekly's to write), and the 11-05 night goes to a priced-family variant or unspent. Trailing-7-day count on 11-05 (cutoff 10-30): 11-01, 02, 03, 04 = four, five with it.

**Why:** C5 proper (entry 50's rule) is the degree-2 root grid, then the column set. Entry 53 brackets the root at 1.25 pF; this window reads it at entry 51's resolution, so C19 (C6 and I36 at degree 2) and `ANS` work at the production order sit on a root as well-read as the degree-1 one. The bracketing pair's interpolation is the reading (entry 51's finding).

**Price (measured in the family):** entry 50 on the same mesh, module, degree and width — build 341.8 s, Z50 pair 2 892.1 s, node pair 3 056.9 s, `memory.peak` 472.42 GiB; entry 53 (three nodes) re-measures the per-node figure on 11-04 and this entry inherits it. **Predicted ≈ 12 400 s, 472–490 GiB; readout brackets window [9 600, 14 400] s, `memory.peak` [440, 512] GiB** (entry 53's). Timeout 14 400 s; 512 GiB limit; the OOM clause as entry 53's.

**Readout (record; C4a's gates asserted at their imported bands as in entries 50 / 53 — (A0), census 64, U₂ 7 885 674 exact, (Deg), kinds and (P_RC) at each node, (P50) and NEG on the Z50 pair; (Z50) printed only):** Im Y_in at each node on all four series; the crossing on P33 (1, +1j) by the bracketing pair, beside entry 53's coarse interpolation and the degree-1 root 193.5996 pF (the order move in pF and %); the degree-2 slope at the root; per-node times and `memory.peak`.

**Decision rule (pre-registered; the first weekly after the window reads it):** every gate green and one crossing inside the grid ⇒ **C_in(F-clinical-1.5T, d2) exists at 0.25 pF resolution** (a record, an internal tuning statement at the production order; no absolute or external claim) and C19 is posed at it; every gate green and no crossing ⇒ entry 53's bracket was off by more than 0.25 pF — recorded, one re-centre, no second window before the weekly; a gate red ⇒ a known-issues entry on the degree-2 fixture; killed / OOM ⇒ the node JSONs written stand as records, the step returns to the weekly. No band moves.

**Measured price (written 2026-11-04 10:30 review from entry 53, the same mesh, module, degree and `-n 16`):** build 368.2 s, Z50 pair 3 149.8 s, node pairs 2 970.4 / 3 258.1 / 3 582.0 s, window 13 351 s, `memory.peak` 480.47 GiB (`20261104T080008Z_GEO-38-C5b.log:13041`, `:13095–13096`). Two nodes ⇒ **predicted ≈ 9 460–10 690 s, 472–490 GiB; readout brackets window [8 000, 13 000] s, `memory.peak` [440, 512] GiB.** Timeout 14 400 s; 512 GiB limit (entry 53 used 93.8 % of it); the OOM clause as entry 53's.

**Command** (entry 53's durable-capture shape, only the node list and the raw-log name changed; the clerk copied it verbatim):

```
XL_CHUNK="GEO-38-C5b-prime"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c5b-prime-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C4=1 FEM_EM_GEO38_C4_NODES_PF=185.6250,185.8750 FEM_EM_GEO38_C4_DEGREE=2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs are gitignored.

### 57. `xl` (Sunday 2026-11-08) — `GEO-38` C16r: entry 49's root bracket refined to 0.05 pF on F-clinical-3T, degree 1, `-n 16` (daily 2026-11-01 03:00, step 6b rule 4 — **priced-family variant** of entry 49, one named knob varied: `FEM_EM_GEO38_C12_NODES_PF`)

**Status:** **READ 2026-11-09 (weekly) — the first branch:** every gate green ((A0), (C), (W), (Az), (U1), (F), (Deg), (P50) 3.921e-13 with NEG 110.1×, (Z50) 6.2e-14, (P_RC) ≤ 2.9e-12 at the four nodes) and one − → + crossing of Im Y_in on P33 (1, +1j) inside the grid — −6.661459e-05 / +1.272741e-03 S at 44.0500 / 44.1000 pF, **44.052487 pF by the bracketing pair** (`:13229`), 0.009 pF (0.021 %) below entry 49's chord, the direction the convexity predicted. All four series cross inside the grid (P41 (1, +1j) 44.078137, P33 (1, −1j) 44.045712, P41 (1, −1j) 44.079569 pF; spread 0.0339 pF against entry 49's 0.0391), a record. ⇒ **C_in(F-clinical-3T, d1) = 44.0525 pF at 0.05 pF resolution** (four decimals, entry 51's convention) — an internal tuning statement at degree 1, no absolute or external claim; the node entries 58 / 60 / 61 carry. Price inside both brackets (1 121 s, 38.11 GiB). No band moves. *Was:* **RUN `20261108T080008Z_GEO-38-C16r.log`** (2026-11-08 02:00; Status 0, 1 121 s, 38.11 GiB; every gate green; one − → + crossing on P33 (1, +1j) at 44.052487 pF by the bracketing pair, `:13229` — the decision rule's first branch as the 03:00 review reads it from a `log-pathologist` ruling, COUNTABLE; ledger row filled; the node is written as 44.0525 pF for entries 58 and 61; the 11-09 weekly reads it). **KEPT 2026-11-02 (weekly): entry 49 read as its first branch; the grid spans all four series' chords and the quadratic estimate, and entry 58 (the 3 T C6) takes this window's node.** **QUEUED 2026-11-01 (03:00 review) for Sunday 2026-11-08 02:00** — `docs/testing/xl-queue.d/20261108-GEO-38-C16r.env` (the command below, copied verbatim with the Write tool). Both lines of the file match this entry. Node forms checked first: `20261101T090538Z_LEAK-CHECK-C16R.log`, clean, Status 0. Trailing-7-day count on 11-08 (cutoff 11-02): 11-02, 03, 04, 05, 06 = five, six with it. The launcher takes the lexically first file, so if the 11-05 or 11-06 night has no file of its own, this window runs that night instead. **The weekly may withdraw it by deleting the file.**

**Why:** the first open `xl` night, and the floor stood at three ahead. Entry 49 found the 3 T root inside a 0.2781 pF bracket, 0.63 % of the root, where entry 51 read the 1.5 T root at 0.13 %; and the weekly's finding on that entry (a bracket's chord is off by the curve's convexity, 0.042 % at 1.5 T over a 0.78 % bracket) applies here: Im Y_in's segment slopes fall 3.54e-2 → 2.59e-2 → 1.83e-2 S/pF across the crossing, so a chord over the bracket reads high. A quadratic through entry 49's three nearest nodes on P33 (1, +1j) reads ≈ 44.05 pF against the chord's 44.0615 (review arithmetic, not a measurement). A 3 T C6 window (after C20) needs one node, and this is the same refinement entry 51 gave the 1.5 T twin before its windows were written.

**Price (measured, entry 49, same module, mesh and width):** build 327.6 s, Z50 pair 135.1 s, node pairs 146.3–173.0 s, `memory.peak` 39.39 GiB. Four nodes ⇒ **predicted ≈ 1 100–1 180 s, 39–42 GiB; readout brackets window [800, 2 700] s, `memory.peak` [30, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Grid (4 nodes, written by this review):** 44.0000, 44.0500, 44.1000, 44.1500 pF, 0.05 pF spacing (0.11 %), none equal to an entry-49 node, spanning all four series' chord readings 44.0505–44.0897 and the quadratic estimate. Node JSONs do not collide.

**Readout (record; asserted at their imported bands as in entry 49 — (A0), census 64, U₁ 1 425 765, (Deg), (Z50) and (P50) and NEG on the Z50 pair, (P_RC) at each node):** Im Y_in at each node on all four series; the crossing on P33 (1, +1j) by the bracketing pair beside entry 49's 44.061525 pF and C15's 43.9719 pF; the four series' spread; per-node times and `memory.peak`.

**Decision rule (pre-registered; the first weekly after the window reads it):** every gate green and one − → + crossing on P33 (1, +1j) inside the grid ⇒ **C_in(F-clinical-3T, d1) at 0.05 pF resolution is the record** (an internal tuning statement, no absolute or external claim) and the weekly poses the 3 T C6 window at it; every gate green and no crossing inside the grid ⇒ entry 49's bracket stands as the record, the miss recorded, no second window; a gate red ⇒ a known-issues entry on the clinical-scale ESR sheet at 128 MHz. No band moves.

**Command** (entry 49's durable-capture shape, only the node list and the raw-log name changed):

```
XL_CHUNK="GEO-38-C16r"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c16r-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C12=1 FEM_EM_GEO38_C12_NODES_PF=44.0000,44.0500,44.1000,44.1500 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_3t_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs are gitignored.

### 58. `xl` (Monday 2026-11-09 at the earliest) — `GEO-38` C6-3T: the ASTM-box-loaded F-clinical-3T under the two-port drive at the degree-1 root, no rod — `|B₁⁺|` and the located 10 g hotspot at 128 MHz, `-n 8` (weekly 2026-11-02, §10 chain C, entry 58; entry 52's family on the 3 T twin — C20's mesh with the rod mapped to gel)

**Status:** **READ 2026-11-09 (weekly) — the first branch: every assertion green ⇒ C6-3T lands as a record.** Status 0, 780 s, `memory.peak` 33.48 GiB at `-n 8` (`:12260`: build 448.3 s, columns 194.1 + 13.4 s, readouts 117.1 s); (A0) 1 128 430 cells rel 0, (U1) 1 508 256 exact, (A0m) census (tag 403 = 14 848 cells, gel-mapped), (R1) 1.893e-13 against 1e-3, (P) 5.120e-13 / 1.810e-12 with NEG 1.2375e-02 / 1.2827e-02 on P33 / P41, (P_RC) 5.120e-13, coverage (i) power 2.975e-14 / mass 3.220e-15 (`:12227–12257`). **Records (degree 1, at the 3 T root 44.0525 pF, V_th = 1 V per feed):** mean `|B₁⁺|`_ccw **1.174932617e-08 T over the 18-point lattice** (the registered readout) and 1.562098621e-08 T over 400 tag-3 centroids, both printed (`:12235`); mean `|B₁⁻|` 9.328405808e-09 T, (Sense) 1.6746, CV 0.3593 (`:12237`); **c\* = x+105 mm dz+50 mm (r 0.105, dz +0.050, az 0°), SAR10g 1.748477085e-04 W/kg (ccw), the cw ball 1.427632775e-04** (`:12256`), the lattice argmax (the box-wide maximum is C6b's reader, off in this window); P_box 3.659141247e-03 / P_src 5.347172304e-02 / P_acc 1.874513035e-02 W; per accepted watt SAR10g(c\*) 9.327634e-03 W/kg/W and centroid `|B₁⁺|`/√P_acc 1.140944e-07 T/√W (`:12259`); S₂ (`:12227`); Im Y_in at the node printed, both senses (`:12233–12234`; ccw P33 +1.032e-03, P41 −1.347e-03 S — opposite signs, a record of the box's detuning, never a red). **Beside entry 52 (C33's reader, `20261109T103557Z_GEO-38-C33-W1-r2.log:194–203`), 3 T / 1.5 T per accepted watt:** lattice mean `|B₁⁺|`/√P_acc **1.410682**, centroid mean 1.344035, **SAR10g(c\*) 4.101440**, P_box / P_acc **3.693500**; c\* moves from x−105 mm dz−50 mm to x+105 mm dz+50 mm — the mirror ball of a near-tie (the top two 1.5 T balls differ by ≈ 6.3e-04 relative), not a moved hot region (at the fixed 1.5 T ball the ratio is ≈ 4.02). **The SAR lever is mostly the box's power share (3.69 of the 4.10).** An internal statement at degree 1: no tuning, compliance or homogeneity claim, no external check. **Posed tonight, by this rule: the 3 T `ANS` case — `ANS-12` (SPEC v1 skeleton; the node TBD = the 3 T degree-2 root, entry 64 then its refinement; our adjudication column the 3 T degree-2 C6 / I36 set, entry 66).** Price inside both brackets; `memory.peak` 33.48 GiB sits 6.5 GiB under the 40–50 GiB prediction (the 3 T box mesh is 7 % smaller and the degree-1 column factorises in 194 s), a price reading. C33 / C47 proceed; entry 60 is READY at this node (below). *Was:* **RUN `20261109T080008Z_GEO-38-C6-3T.log`** (2026-11-09 02:00). Status 0 at 780 s, `memory.peak` 33.48 GiB. Every gate is green: 5 passed and 7 skipped on all 8 ranks by token, orphans 0 / 0. The 11-09 weekly reads the decision rule. The 2026-11-09 03:00 review's mechanical read is in the ledger row. *Was:* **QUEUED 2026-11-08 (03:00 review) for Monday 2026-11-09 02:00** — `docs/testing/xl-queue.d/20261109-GEO-38-C6-3T.env`, the command below with `<node>` = **44.0525 pF** (entry 57's crossing 44.052487 pF on P33 (1, +1j), `20261108T080008Z_GEO-38-C16r.log:13229`, at four decimals), written with the Write tool and diffed identical to the block below after the substitution. Node forms (44.0525, its `.12e` pF and F forms, the five-place crossing) cleared by the leak guard: `20261108T090440Z_LEAK-CHECK-C6-3T.log`, Status 0. Trailing-7-day count on 11-09 (cutoff 11-03): 11-03, 04, 05, 06, 08 = five, six with it. **READY EXCEPT ITS NODE 2026-11-04 (weekly): (i) C21 MET, (iii) entry 52 read as its first branch tonight; only (ii) entry 57 (Sunday 11-08) is open.** The first review after entry 57 substitutes `<node>` below with 57's crossing on P33 (1, +1j) by the bracketing pair (or 44.0615 pF if 57 is withdrawn or finds no crossing inside its grid), checks the node's printed forms against the leak guard, and queues `xl-queue.d/20261109-GEO-38-C6-3T.env` for Monday 11-09 02:00 (trailing-7-day count on 11-09, cutoff 11-03: 11-03, 04, 05, 06, 08 = five, six with it). **Price ratified from entries 49 and 52 (the post-processing half now measured):** build ≈ 460 s (C20 433 s of gmsh; C21's build window 472 s at `-n 2`), one factorising column with two held ≈ 250 s (entry 48's arm 263 s on 1 217 025 cells at `-n 8`; entry 49's node pair 146–173 s at `-n 16`), readouts ≈ 170 s (entry 52's 169.7 s) ⇒ **predicted ≈ 890–1 000 s, 40–50 GiB; readout brackets window [700, 1 800] s, `memory.peak` [35, 70] GiB.** **Readout as entry 52's with its two corrections:** the 18-point lattice mean and the 400-centroid mean are named separately (the lattice mean is the registered readout); c\* is the lattice argmax, not a box maximum (C6b's box-wide reader is not in this window). **Command** (entry 52's durable-capture shape with C21's module and knob; the clerk substitutes `<node>` and copies it verbatim):

```
XL_CHUNK="GEO-38-C6-3T"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c6-3t-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C21=1 FEM_EM_GEO36_S8_DRIVE_C_PF=<node> timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_3t_box_loaded.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

*Was:* **PENDING PREREQUISITES (written 2026-11-02 weekly):** (i) **§10 C21** — the 3 T box-loaded body, a slot step: C6a's module (`tests/validation/test_geo38_f_clinical_box_loaded.py`) parametrised on the 3 T twin (C20's mesh, 1 128 430 cells, U₁ 1 508 256, census 64; 128 MHz; the node from a required knob; r = 0.2 Ω on all 64 gaps; the two-port quadrature drive on C16's feeds P33 / P41), build-only green at `-n 2` and the K-gates stub-run (items 701 / 702's pattern). **(i) MET 2026-11-03 (§9 items 775 + 791):** the module `tests/validation/test_geo38_f_clinical_3t_box_loaded.py` (knob `FEM_EM_GEO38_C21`, node knob `FEM_EM_GEO36_S8_DRIVE_C_PF`) stub-ran green (`20261103T103420Z_GEO-38-C21-W0.log`), and its `build` window is green at `-n 2`: (A0) 1 128 430 cells, rel 0, U₁ 1 508 256 exact, census 64, (A0m) none, 18 / 18 lattice points in tag 3, (Bore) gap 0.137232684 m, 5 passed per rank, window 472 s (`20261103T193942Z_GEO-38-C21-W1-build-r2.log`); (ii) **entry 57 green** (Sunday 11-08) — its crossing is the node, or 44.0615 pF (entry 49's) if 57 is withdrawn or finds no crossing inside its grid; (iii) **entry 52 read** by the 11-04 weekly (the 1.5 T body's first unstubbed run — a red there is read before the 3 T twin is bought). **The first review after all three writes the command in entry 52's shape with the 3 T module / knob and the node, `-n 8`, marks READY and queues it into `xl-queue.d/20261109-GEO-38-C6-3T.env`** (trailing-7-day count on 11-09, cutoff 11-03: 11-03, 04, 05, 06, 08 = five, six with it); node forms through the leak guard first. If (i) or (ii) is missing by 11-08 18:00, the 11-09 night goes to a priced-family variant under the daily's licence or unspent.

**Why:** the operator's ladder step (6) on the 3 T twin — the ASTM box loaded, `|B₁⁺|`, 10 g SAR, a located hotspot at 128 MHz — on the fixture the 2026-10-14 directive names second; the first 128 MHz clinical-scale safety quantity, and the fixture a 3 T `ANS` case is posed on. C20's rod-in-box mesh is reused with the rod mapped to gel (I41's route), as entry 52 reuses I36's, so no new mesh is bought.

**Price (measured in the family; scaling stated, the post-processing half to be read from entry 52):** entry 49 on the 3 T coil mesh (1 057 895 cells) at `-n 16`: build 327.6 s, Z50 pair 135.1 s, node pairs 146.3–173.0 s, 39.39 GiB; C20's box mesh is 1 128 430 cells (+6.7 %), one factorisation with two held columns at `-n 8` ≈ entry 48's arm (263 s on 1 217 025 cells at `-n 8`) ⇒ ≈ 250 s; build ≈ 460 s (C20: 463 s at `-n 2`); the lattice and balls as entry 52 measures them (predicted there 82 s). **Predicted ≈ 850–1 000 s, 40–50 GiB; readout brackets window [600, 2 400] s, `memory.peak` [30, 120] GiB** (entry 52's, until 52 is read). Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands as entry 52's — (A0) cells 1 128 430 within `CELL_COUNT_BAND`, (A0m) the matched census exact, census 64, U₁ exact; (P) and NEG; (R1) on the 2×2 at `RECIPROCITY_BAND`; (P_RC) at the node at `DISCRETE_IDENTITY_RTOL`; coverage (i) power / mass at 1e-10):** mean \|B₁⁺\|_ccw over the 18-point lattice and its z_c spread; SAR10g(c\*) and its location; P_box / P_src and the per-box-watt figures beside the per-accepted-watt ones (the normaliser stated); Im Y_in at the node (printed, never a red); per-phase time and RSS. No mirror gate, no class gate: the maps are records.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ C6-3T lands as a record — the first 128 MHz clinical-scale `|B₁⁺|` / 10 g figure at degree 1 at the root, an **internal statement** — printed beside entry 52's 1.5 T figures (the 1.5 T / 3 T ratio of the mean and of the hotspot, a record); the 3 T `ANS` case is posed by that weekly beside C7; (P) / (R1) / NEG red with the identities green ⇒ a known-issues entry on the clinical-scale conductor model at 128 MHz with the box load, the maps as records; (A0) / (A0m) or an identity red ⇒ the step returns to the weekly, no second window; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command:** written by the READY review in entry 52's durable-capture shape (`R=/workspace/logs/geo38-c6-3t-xl-raw.log`, the 3 T module and knob, `FEM_EM_GEO36_S8_DRIVE_C_PF=<node>` or the knob C21 names, `-n 8`). No `XL_OUTPUTS`: the maps and the node JSON are gitignored.

### 59. `xl` (the first open night after entry 58, 2026-11-10 at the earliest) — `GEO-38` P33: the clinical box heated by the coil — C6's SAR heats the ASTM gel box inside F-clinical-1.5T, no implant, degree 1, `-n 8` (weekly 2026-11-02, §10 chain P, entry 59; entry 52's family — the same mesh, the same drive, `MAT-15`'s thermal route added, as entry 55 adds it to entry 48)

**Status:** **READ 2026-11-11 (weekly): the first branch — every *executed* assertion green, so P33 lands as a record (θ(t\*) at c\* per accepted watt 1.611905088e-03 K ccw, `:12326`) and the 1.5 T thermal `ANS` case is posed (`ANS-13`, SPEC v1); the readout's "hotspot-ball slope" and "(S0)" were never in the module — the 11-02 readout wording outran P61 — and the DG0 tip print is known-issues 2026-11-06's scope; `cem43` per watt is the 20 °C baseline, not a dose. No band moves. The degree-2 twin is entry 69.** *Was:* **RUN `20261110T080008Z_GEO-38-P33.log`** (2026-11-10 02:00; Status 0, 852 s, 41.70 GiB at `-n 8`; 5 passed / 9 skipped on all 8 ranks by token; every executed assertion green, EM and thermal; **the readout's hotspot-ball slope and "(S0)" were never asserted under those names** — the module asserts the global (S2) only (known-issues 2026-11-06 scope) — and the hotspot `cem43` per accepted watt sits at the 20 °C baseline; the 2026-11-10 03:00 review's mechanical read from a `log-pathologist` ruling is in the ledger row; **the 11-11 weekly reads the decision rule**, including whether "every assertion" means the executed ones). *Was:* **QUEUED 2026-11-05 (10:30 review) for Tuesday 2026-11-10 02:00** — `docs/testing/xl-queue.d/20261110-GEO-38-P33.env`, the command below copied verbatim with the Write tool (diffed identical); it already carries the landed knob `FEM_EM_GEO38_C6_HEAT_SENSE=both`. Node forms (193.5996, its `.12e` pF and F forms) cleared by the leak guard: `20261105T163300Z_LEAK-CHECK-P33.log`, Status 0. Trailing-7-day count on 11-10 (cutoff 11-04): 11-04, 05, 06, 08, 09 = five, six with it. **Lexical-order note:** if entry 58's `20261109-…` file is not written before the 11-09 02:00 window, the launcher takes this file that night instead (no dependency on entry 58 is lost; the 11-08 03:00 review writes 58's file). **(P61b) MET 2026-11-05 (§9 item 818, landing commit of `test(geo-38): P61b`; `20261105T103249Z_GEO-38-P61b-W0.log`, `20261105T103333Z_GEO-38-P61b-W0-n1.log`; knob `FEM_EM_GEO38_C6_HEAT_SENSE` = `both` / `ccw` / `cw`, default `cw` — the landed name is the one in the command below).** **READY ON P61b 2026-11-04 (weekly): (i) entry 52 read as its first branch tonight, (ii) P61 MET; the sense is ruled *both*.** P61 heats the (1, +1j) sense (cw, entry 55's drive) at C6's c\*, which was located under ccw (1, −1j); the thermal arm costs seconds, so **P61b** (§10 chain P, a no-solve slot step) adds `FEM_EM_GEO38_C6_HEAT_SENSE` = `both` / `ccw` / `cw` and runs both arms in one window, each per its own P_acc, the ccw figures (the record's sense) reported first and cw beside them. **The first review after P61b lands writes the command below with the knob, checks the node's printed forms against the leak guard and queues `xl-queue.d/20261110-GEO-38-P33.env` for Tuesday 11-10 02:00** (trailing-7-day count on 11-10, cutoff 11-04: 11-04, 05, 06, 08, 09 = five, six with it). **Price from entry 52 (measured) plus two thermal arms:** EM half 868 s, 40.95 GiB (`20261103T080009Z_GEO-38-C6.log:12312`, `:12404`); each thermal arm one CN transient to t\* on ≈ 100 000 gel cells plus the source, ball and probe reads, ≈ 30–120 s (entry 55's pricing from I41 and entry 42) ⇒ **predicted ≈ 950–1 150 s, 41–50 GiB; readout brackets window [800, 1 800] s, `memory.peak` [35, 70] GiB.** Timeout 14 400 s; 512 GiB limit. **Command** (entry 52's durable-capture shape with P61's module and knobs; the clerk copies it verbatim once P61b lands — if P61b's knob name differs, the clerk uses the landed name and says so):

```
XL_CHUNK="GEO-38-P33"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-p33-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C6_HEAT=1 FEM_EM_GEO38_C6_HEAT_SENSE=both FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_box_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

*Was:* **PENDING PREREQUISITES (written 2026-11-02 weekly):** (i) **entry 52 green** (Tuesday 11-03; its first branch: every assertion green ⇒ C6 lands as a record and P33 is posed); (ii) **§10 P61** — the heating body, a slot step: `tests/validation/test_geo38_f_clinical_box_heating.py` (knob `FEM_EM_GEO38_C6_HEAT=1`, the node knob as C6's), C6's field under the (1, +1j) drive at the node, its σ|E|²/ρ on the gel (tags 3, 403) as the DG0 source of `MAT-15`'s thermal route on the same mesh, Robin h = `MAT-21`'s on the six faces, CN to t\* with step 1's Δt rule; stub-run on a small `create_box` with `_drive` stubbed (I32a's pattern, `20261101T103929Z_GEO-36-I32a-W0b.log`). **(ii) MET 2026-11-03 (§9 item 777):** the module landed with every K-gate green and both named controls firing, 7 passed / 1 skipped per rank at `-n 2` and `-n 1` (`20261103T104609Z_GEO-38-P61-W0c.log:138`, `20261103T104625Z_GEO-38-P61-W0c-n1.log:117`). It heats the (1, +1j) sense with the ccw P_gel printed beside it, and reports θ, `cem43` and the 360 s running mean per accepted watt (P_acc). Which sense entry 59 reports is the weekly's call. **The first review after both land (the 11-04 weekly reads entry 52) writes the command in entry 52's shape with `FEM_EM_GEO38_C6_HEAT=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996`, `-n 8`, prices it from entry 52's measured window plus the thermal phase (seconds, as entry 55's pricing from I41 and entry 42), marks READY and queues it into `xl-queue.d/20261110-GEO-38-P33.env`** (or the first open night after entry 58; trailing-7-day count on 11-10, cutoff 11-04: 11-04, 05, 06, 08, 09 = five, six with it). If entry 52 is red or killed, this entry waits for the weekly's re-pose; if P61 has not landed by 11-09 18:00, the night goes to a priced-family variant or unspent.

**Why:** subgoal (8e)'s FEM half — the ASTM F2182 gel box heated by the clinical coil, by the route `MAT-15` gated on F-small and `ANS-9` checked externally (one-way, gel): the first clinical-scale heating figure without an implant, beside entry 55's with one. The closed-form half is `MAT-21` (P1 … P60); per accepted watt on a 20 °C gel; a record: no tissue, no compliance claim, no external check.

**Price (to be filled from entry 52; the EM half measured, the thermal half scaled):** the EM half is entry 52's window (predicted ≈ 850–1 000 s, 42–50 GiB at `-n 8`); the thermal half is one CN transient to t\* on the gel (P1, DG0 source) — I41's 2.5 / 1.7 s on 66 693 / 51 693 sub cells and entry 42's 0.55–0.84 s at `-n 16` are the measured points, ≈ 30–120 s with the source, ball and probe reads. **Readout brackets, provisional until priced:** window [900, 2 400] s, `memory.peak` [35, 120] GiB. Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands — C6's EM gates ((A0), (A0m), census 64, U₁, (P) / NEG, (R1), (P_RC), coverage (i)); the scheme's discrete energy identity ≤ `E_BAND`; the calorimetric initial slope dT/dt|₀ = SAR/c at `MAT-15`'s band on the gel and on the 10 g hotspot ball; (S0) census):** θ at the hotspot ball centre and at P37's six box probe points at t\*; the six-minute window and `cem43` at the hotspot (`MAT-20`'s forms, per accepted watt, the normaliser stated); per-phase time and RSS.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ P33 lands as a record (no band on θ) and the 1.5 T `ANS` thermal case on F-clinical (HFSS + Icepak, `ANS-9`'s pattern at clinical scale) is posed by that weekly beside C7; an energy or slope identity red with the EM gates green ⇒ a known-issues entry on the thermal route at clinical scale, the EM record stands; an EM gate red ⇒ as entry 52's rule; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command:** written by the READY review in entry 52's durable-capture shape with P61's module and knobs (`R=/workspace/logs/geo38-p33-xl-raw.log`, `-n 8`). No `XL_OUTPUTS`: the JSON and the fields are gitignored.

### 60. `xl` (Wednesday 2026-11-11 at the earliest) — `GEO-38` P33-3T: the 3 T clinical box heated by the coil — C6-3T's SAR heats the ASTM gel box inside F-clinical-3T at 128 MHz, no implant, degree 1, `-n 8` (weekly 2026-11-04, §10 chain C / P, entry 60; entry 59's family on the 3 T twin — entry 58's mesh, drive and node, `MAT-15`'s thermal route added)

**Status:** **READ 2026-11-11 (weekly): the first branch — every executed assertion green, so P33-3T lands as a record (θ(t\*) at c\* per accepted watt 6.516151132e-03 K ccw, `:12276`; the 3 T / 1.5 T lever 4.0425 ccw / 4.5397 cw by C39) and the 3 T thermal case is posed beside the 1.5 T one (`ANS-13`'s second column, its node `ANS-12` v2's); the gel's constants are the 1.5 T values (provisional, W3). No band moves. The degree-2 twin is entry 70.** *Was:* **RUN `20261111T080008Z_GEO-38-P33-3T.log`** (2026-11-11 02:00; Status 0, 732 s, 32.75 GiB at `-n 8` — under the [35, 70] GiB readout bracket, a price reading; 5 passed / 5 skipped on all 8 ranks by token; every executed assertion green, EM and thermal; **as on entry 59, the readout's hotspot-ball slope and "(S0)" were never asserted under those names**, and the hotspot `cem43` per accepted watt sits at the 20 °C baseline; the 2026-11-11 03:00 review's mechanical read from a `log-pathologist` ruling is in the ledger row; **the 11-11 weekly reads the decision rule**). C39 (§9 item 953) reads it beside entry 59. *Was:* **READY and QUEUED 2026-11-09 (weekly) for Wednesday 2026-11-11 02:00** — (i) C31 MET (below), (ii) entry 58 read tonight on its first branch (above). `docs/testing/xl-queue.d/20261111-GEO-38-P33-3T.env`, the command below, written with the Write tool: entry 59's command with the raw-log name, the module (C31's `tests/validation/test_geo38_f_clinical_3t_box_heating.py`), the knob (`FEM_EM_GEO38_C31=1` for `FEM_EM_GEO38_C6_HEAT=1`) and the node (44.0525) changed, `FEM_EM_GEO38_C6_HEAT_SENSE=both` kept (P61b's knob, which C31 shares). Node forms: 44.0525 and its printed forms were cleared for entry 58 (`20261108T090440Z_LEAK-CHECK-C6-3T.log`, Status 0) — the same node, no new form. Trailing-7-day count on 11-11 (cutoff 11-05): 11-05, 06, 08, 09, 10 = five, six with it. **Lexical-order note:** this file sorts before entry 61's `20261112-…`, so 11-11 is this window's and 61 keeps 11-12. **The clerk diffs the file against this entry once.** Price unchanged from the pre-registration (predicted ≈ 950–1 200 s, 40–50 GiB; brackets [700, 1 800] s, [35, 70] GiB) — entry 58 measured 780 s / 33.48 GiB for the EM half, so the thermal arms have ≈ 1 000 s of bracket.

```
XL_CHUNK="GEO-38-P33-3T"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-p33-3t-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C31=1 FEM_EM_GEO38_C6_HEAT_SENSE=both FEM_EM_GEO36_S8_DRIVE_C_PF=44.0525 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_3t_box_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

*Was:* **PENDING PREREQUISITES (written 2026-11-04 weekly):** (i) **§10 C31** — P33-3T's body, a slot step: P61's module parametrised on C21's 3 T twin (C21's build, 128 MHz, the node knob, both senses as P61b's), stub-run with P61's K-gates and both controls — **(i) MET 2026-11-05 (§9 item 842; module `tests/validation/test_geo38_f_clinical_3t_box_heating.py`, gated test `test_geo38_c31_box_heating_3t`; knobs `FEM_EM_GEO38_C31=1` plus the node `FEM_EM_GEO36_S8_DRIVE_C_PF` (required) and the sense `FEM_EM_GEO38_C6_HEAT_SENSE` (P61b's); degree 1 only — `FEM_EM_GEO36_S8_DRIVE_DEGREE` is refused by its foreign guard; `20261105T181050Z_GEO-38-C31-W0r.log`, `20261105T181131Z_GEO-38-C31-W0r-n1.log`)**; (ii) **entry 58 read** by the first weekly after it (Wednesday 11-11 at the earliest). **The first review after both writes the command in entry 59's shape with C31's module and knobs and entry 58's node, marks READY and queues it into `xl-queue.d/20261111-GEO-38-P33-3T.env`** (trailing-7-day count on 11-11, cutoff 11-05: 11-05, 06, 08, 09, 10 = five, six with it); node forms through the leak guard first. If either is missing by 11-10 18:00, the night goes to a priced-family variant under the daily's licence or unspent.

**Why:** subgoal (8e) on the fixture the 2026-10-14 directive names second — the first 128 MHz clinical-scale heating figure, beside entry 59's at 64 MHz; the 1.5 T / 3 T ratio of the hotspot temperature per accepted watt is the reading. A record: no tissue, no compliance claim, no external check.

**Price (measured in the family):** entry 58's window (predicted ≈ 890–1 000 s, 40–50 GiB) plus two thermal arms of ≈ 30–120 s each ⇒ **predicted ≈ 950–1 200 s, 40–50 GiB; readout brackets window [700, 1 800] s, `memory.peak` [35, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands — entry 58's EM gates; the scheme's discrete energy identity ≤ `E_BAND`; the calorimetric initial slope at `MAT-15`'s band on the gel and on the 10 g hotspot ball; (S0) census):** θ at the hotspot ball centre and at P37's six box probe points at t\* under both senses; the six-minute window and `cem43` at the hotspot per accepted watt (the normaliser stated); the 1.5 T / 3 T ratio beside entry 59's figures; per-phase time and RSS.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ P33-3T lands as a record (no band on θ) and the 3 T thermal `ANS` case is posed beside the 3 T EM one; an energy or slope identity red with the EM gates green ⇒ a known-issues entry on the thermal route at 128 MHz, the EM record stands; an EM gate red ⇒ as entry 58's rule; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command:** written by the READY review in entry 59's shape (`R=/workspace/logs/geo38-p33-3t-xl-raw.log`, C31's module and knobs, `-n 8`). No `XL_OUTPUTS`.

### 61. `xl` (Thursday 2026-11-12 at the earliest) — `GEO-36` step 11 / I36-3T: the F2182 rod in the ASTM box inside F-clinical-3T under the two-port drive at the degree-1 root, 128 MHz, `-n 8` (weekly 2026-11-04, §10 chain C / I, entry 61; entry 48's family on the 3 T twin — C20's rod-in-box mesh, tag 403 kept as σ_wire cells)

**Status:** **READ 2026-11-14 (weekly): the first branch — every asserted gate green ⇒ I36-3T lands as a record (tip 1 g Q 4.748393 / 4.884059 top / bottom, 4.748393 / 4.775105 apices; η_rod(T1) 3.407173; rod share of the volume power 7.434370e-05; Im Y_in P33 +1.880e-03 / P41 −1.974e-03 S, signs opposite to entry 48's — a record) at the 3 T degree-1 root; against entry 48 by C40 (§9 item 967): Q per ball ≈ 1.08–1.10×, η_rod(T1) 0.935×, rod share 1.99× — the frequency lever on one rod model and one definition, at the gel's 1.5 T σ / ε_r (provisional, W3). By the rule the 3 T implant-heating window is posed — entry 67, pre-registered 11-09, re-dated tonight behind entry 71. The 34.23 GiB under the bracket is a price reading. No band moves.** *Was:* **RUN 2026-11-12 `20261112T080008Z_GEO-36-step11-I36-3T.log`** — Status 0, 848 s, `memory.peak` 34.23 GiB (under the [35, 70] GiB bracket, a price reading), 5 passed / 4 skipped on all 8 ranks, every asserted gate green (ledger row, `log-pathologist` reading, 2026-11-12 03:00 review). Readout records: tip 1 g Q 4.748393–4.884059, η_rod(T1) 3.407173, rod share 7.434370e-05. **The 11-14 weekly reads the decision rule;** the mechanical read opens C40 (§9 item 967) under the after-window licence and makes entry 67's (ii) MET. *Was:* **READY and QUEUED 2026-11-08 (03:00 review) for Thursday 2026-11-12 02:00** — (i) C32 MET (below), (ii) entry 57 green (`20261108T080008Z_GEO-38-C16r.log`, the first branch; node 44.0525 pF). `docs/testing/xl-queue.d/20261112-GEO-36-step11-I36-3T.env`, the command below, written with the Write tool and diffed against entry 48's command with only the raw-log name, the module (C32's), the knob (`FEM_EM_GEO36_S11_DRIVE=1` for `FEM_EM_GEO36_S8_DRIVE=1`) and the node changed. Node forms cleared: `20261108T090440Z_LEAK-CHECK-C6-3T.log`. Trailing-7-day count on 11-12 (cutoff 11-06): 11-06, 08, 09, 10 = four, plus 11-11 if a window runs then — six at most with it. **Lexical-order note:** the 11-11 night has no file (entry 60 waits on the 11-11 weekly's reading of entry 58), so the launcher takes this file on 11-11 unless a review writes a `20261111-…` file first; nothing is lost either way (entry 61 waits only on 57 and C32).

```
XL_CHUNK="GEO-36-step11-I36-3T"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo36 && R=/workspace/logs/geo36-i36-3t-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S11_DRIVE=1 FEM_EM_GEO36_S8_DRIVE_C_PF=44.0525 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_3t_drive.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

*Was:* **PENDING PREREQUISITES (written 2026-11-04 weekly):** (i) **§10 C32** — I36-3T's body, a slot step: I36b's arm M / arm R route on C20's 3 T rod mesh (1 128 430 cells, U₁ 1 508 256, census 64; 128 MHz; the node knob; the two-port drive on C16's feeds P33 / P41), stub-run as I36c / I36d were — **(C32) MET 2026-11-05 (§9 item 822; module `tests/validation/test_geo36_rod_f_clinical_3t_drive.py`, gated test `test_geo36_c32_rod_drive_3t`; knobs `FEM_EM_GEO36_S11_DRIVE=1` plus the node `FEM_EM_GEO36_S8_DRIVE_C_PF` (required); degree 1 only — `FEM_EM_GEO36_S8_DRIVE_DEGREE` and `FEM_EM_GEO36_S8_DRIVE` are refused by its foreign guard; `20261105T133319Z_GEO-36-C32-W0.log`, `20261105T133349Z_GEO-36-C32-W0-n1.log`)**; (ii) **entry 57 green** (Sunday 11-08) — its crossing is the node, or 44.0615 pF if 57 is withdrawn or finds none. **The first review after both writes the command in entry 48's shape with C32's module and knobs and the node, marks READY and queues it into `xl-queue.d/20261112-GEO-36-step11-I36-3T.env`** (trailing-7-day count on 11-12, cutoff 11-06: 11-06, 08, 09, 10, 11 = five, six with it); node forms through the leak guard first. If either is missing by 11-11 18:00, the night goes to a priced-family variant or unspent.

**Why:** the rod in the clinical box at 128 MHz — the second half of (7b) at clinical scale on the 3 T twin; the tip 1 g Q and η_rod(T1) by entry 48's definition beside entry 48's 4.34–4.46 / 3.642 at 64 MHz (the frequency lever on one rod model and one definition), the rod's power share, Im Y_in at the node. A record: no published comparand until (7c).

**Price (measured in the family):** entry 48 on the 1.5 T rod mesh at `-n 8`: build 461.7 s, arm M 263.5 s, arm R 262.3 s, 42.39 GiB; C20's 3 T mesh is 7 % smaller and its gmsh read 433 s; entry 49's 3 T coil factorisation 146–173 s at `-n 16` ⇒ **predicted ≈ 1 000–1 150 s, 40–50 GiB; readout brackets window [800, 1 800] s, `memory.peak` [35, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands as entry 48's — (A0), (A0m), census 64, U₁, (D-rod) G1 / G2 / G4, (M2), (Z), (P) / NEG on both arms, (R1), (P_RC)):** tip 1 g Q on the four balls and η_rod(T1) beside entry 48's and I57's; the rod's share of the volume power; Im Y_in at the node on both feeds; arm M's S beside arm R's; per-phase time and RSS. Box containment not asserted on this route (entry 48's correction).

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ I36-3T lands as a record and the 3 T implant-heating window (entry 55's twin) is posed by that weekly; (P) / (R1) / NEG red with the identities green ⇒ a known-issues entry on the rod in the 128 MHz box, the maps as records; (A0) / (A0m) or an identity red ⇒ the step returns to the weekly, no second window; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command:** written by the READY review in entry 48's shape (`R=/workspace/logs/geo36-i36-3t-xl-raw.log`, C32's module and knobs, `-n 8`). No `XL_OUTPUTS`.

### 62. `xxl` (Saturday 2026-11-14) — `GEO-38` C19: the degree-2 column set at the degree-2 root — C6 (the box-loaded 1.5 T coil, `|B₁⁺|` / 10 g with the box-wide reader) and I36 (the rod in the box) at degree 2, 64 MHz, `-n 16` (weekly 2026-11-04, §10 chain C19; entries 52 / 48's modules at the order of entries 50 / 53 — the column `ANS-11` adjudicates on)

**Status:** **READ 2026-11-14 (weekly): the green branch — every asserted gate green on both calls ⇒ the degree-2 clinical column exists. `ANS-11`'s `COMPARISON.md` degree-2 rows are filled from the tracked log (C37, §9 item 983, `20261114T181245Z_ANS-11-C37-page-r.log`; (J62) 53 fields, 0 mismatches against the host JSONs — custody only) and the degree-1 → degree-2 moves are records: the per-watt figures move −0.46 / −0.74 / −0.49 / −1.26 % with the node change (193.5996 → 185.7629 pF) included — under the 10-24 exit's 2 % rule, so the rule does not fire: entry 52's degree-1 figures stand as records with the move beside them, and the adjudication column is the degree-2 one as SPEC v2 says. The box-wide 10 g maximum g\* (x+157.5 mm, 2.340 × c\*'s) is a record of the C6b reader at degree 2. The rod's degree-2 reading (call 2) is §10 C56's, against entry 48 under the same rule (predicted to fire there). The measured window (call 1 4 301.9 s, call 2 6 590.9 s, 486.57 GiB) re-prices entry 66 and prices entry 69 tonight. No band moves.** *Was:* **RUN 2026-11-14 `20261114T080009Z_GEO-38-C19.log`** (`3185acfd`; Status 0, `rc1=0` / `rc2=0`, 10 902 s, `memory.peak` 486.57 GiB — inside both brackets; every asserted gate green on both calls; ledger row filled by the 2026-11-14 10:30 review from a `log-pathologist` reading; the decision rule's green branch reads mechanically, so C37 (§9 item 983) opens under the after-window licence; the 11-16 weekly reads the record). *Was:* **READY 2026-11-07 (weekly): (i) entry 56 read — the first branch, the root 185.762904 pF by the bracketing pair on P33 (1, +1j), written as the node 185.7629 pF; (ii) C34 MET; (iii) C6b MET.** The command below is entry 53's durable-capture shape as **two pytest calls in one window** (K-62: C6a's foreign guard refuses I36's knobs — call 1 is C6a's module at degree 2 with the box-wide grid, call 2 is I36b's module at degree 2; each under `timeout -k 60 14400`, the window's 28 800 s ceiling; the window's `rc` is the first non-zero of the two and both are printed), `-n 16`, the node's printed forms cleared by the leak guard (`20261108T030938Z_LEAK-CHECK-C19.log`, clean). **QUEUED 2026-11-07 (weekly) — `xxl-queue.d/20261114-GEO-38-C19.env`, written with the Write tool, identical to the block below** (one `xxl` per trailing seven days: 11-14 is seven days after 11-07; the clerk diffs the file against this entry once). Price unchanged (predicted ≈ 12 000–14 000 s, 500–540 GiB; brackets [9 000, 21 000] s, [440, 700] GiB; the grid's ≲ 610 s at degree 1 and its ≈ 72 ball compiles sit inside the bracket — a shared-compile reader was declined). Readout as written, with entry 52's two corrections (the lattice mean and the centroid mean named separately; the box-wide 10 g maximum from C6b's grid beside the lattice argmax c\*).

```
XL_CHUNK="GEO-38-C19"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 /workspace/logs/geo36 && R=/workspace/logs/geo38-c19-xxl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C6=1 FEM_EM_GEO38_C6_DEGREE=2 FEM_EM_GEO38_C6_GRID=1 FEM_EM_GEO36_S8_DRIVE_C_PF=185.7629 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_box_loaded.py -v -s --tb=short >> \$R 2>&1; rc1=\$?; echo \"[C19] call 1 (C6 degree 2, grid) rc1=\$rc1\" >> \$R; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_DRIVE=1 FEM_EM_GEO36_S8_DRIVE_DEGREE=2 FEM_EM_GEO36_S8_DRIVE_C_PF=185.7629 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/validation/test_geo36_rod_f_clinical_drive.py -v -s --tb=short >> \$R 2>&1; rc2=\$?; echo \"[C19] call 2 (I36 degree 2) rc2=\$rc2\" >> \$R; rc=\$rc1; if [ \$rc2 -ne 0 ]; then rc=\$rc2; fi; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

*Was:* **PENDING PREREQUISITES (written 2026-11-04 weekly):** (i) **entry 56 read** by the 11-07 weekly (the degree-2 root at 0.25 pF by the bracketing pair, or entry 53's chord 185.763346 pF if 56 finds no crossing inside its grid); (ii) **§10 C34** — the degree-2 knob on C6a's and I36b's modules, stub-run, (Deg) U₂ = 2E + 2F exact at build — **(C34) MET 2026-11-05 (§9 item 820, landing commit of `test(geo-38): C34`; `20261105T120358Z_GEO-38-C34-W0.log`, `20261105T120420Z_GEO-38-C34-W0b.log`, `20261105T120440Z_GEO-38-C34-W0c.log`, `20261105T120517Z_GEO-38-C34-W0-n1.log`; knobs `FEM_EM_GEO38_C6_DEGREE=2` (new, C6a) and `FEM_EM_GEO36_S8_DRIVE_DEGREE=2` (I36b's, existing); K-62: C6a's foreign guard refuses I36's knobs, so the command is two pytest calls in one window, each with its own environment prefix — not one call as written below)**; (iii) **§10 C6b** (the box-wide 10 g reader) landed, or the window runs without it and C6b's first window is the next — **(C6b) MET 2026-11-05 (§9 item 844, attempt 2 of 4; `20261105T182202Z_GEO-38-C6b-W0r.log`, `20261105T182247Z_GEO-38-C6b-W0br.log`, `20261105T182320Z_GEO-38-C6b-W0cr.log`, `20261105T182346Z_GEO-38-C6b-W0r-n1.log`); knob `FEM_EM_GEO38_C6_GRID=1` (off by default); measured stub ball time median 0.29 s, max 2.06 s, sum 10.79 s over 31 new balls at `-n 2` — a warm-JIT-cache figure on the stub mesh (attempt 1 had compiled those centres), not a cold compile and not a clinical-mesh price — beside the arithmetic grid cost ≲ 610 s at degree 1, `-n 8`; the degree-2 cost is the weekly's)**. **The 11-07 weekly writes the command (both modules in one pytest call in entry 53's shape with the two degree knobs and the node, `-n 16`), marks READY, and the clerk queues it into `xxl-queue.d/20261114-GEO-38-C19.env`** (one `xxl` per trailing seven days: 11-14 is seven days after 11-07); node forms through the leak guard first. If (i) or (ii) is missing by 11-13 18:00, the Saturday goes to a priced-family variant or unspent.

**Why:** the 10-12 epitaph — a benchmark column is recorded at the production convention (degree 2, exact rectangles) before it is adjudicated; `ANS-11`'s SPEC quotes entry 52's degree-1 column as a record and adjudicates on this one. The rod's degree-2 reading (I36 at degree 2) rides the same window because the factorisation is the price and the arms are two of it.

**Price (measured in the family; scaling stated):** entry 53 on the 1.5 T coil mesh (1 147 366 cells) at degree 2, `-n 16`: build 368.2 s, one factorisation per node pair 2 970–3 582 s (MUMPS 2 596–3 200 s), 480.47 GiB; the box mesh is 1 217 025 cells (+6 %). Three factorisations (C6's column, I36's arm M and arm R) ≈ 3 × 3 300–3 900 s, two builds ≈ 900 s, degree-2 readouts ≈ 3 × 350 s ⇒ **predicted ≈ 12 000–14 000 s, 500–540 GiB; readout brackets window [9 000, 21 000] s, `memory.peak` [440, 700] GiB.** Timeout 28 800 s; 754 GiB limit; an OOM moves the two modules to separate `xxl` Saturdays (the weekly's to write).

**Readout (record; asserted at their imported bands — entry 52's gates for C6 and entry 48's for I36, (Deg) exact; (P_RC) at the node):** entry 52's rows at degree 2 (the lattice mean and the centroid mean separately, c\* and the box-wide 10 g maximum and location, P_box / P_acc, S and Y_in at the feeds) and entry 48's (tip 1 g Q, η_rod(T1), the rod's power share); the degree-1 → degree-2 move of each beside it; per-phase time and `memory.peak`.

**Decision rule (pre-registered; the 11-16 weekly reads it):** every assertion green ⇒ the degree-2 column exists — `ANS-11`'s `COMPARISON.md` rows are filled from it, the degree-1 → degree-2 moves are records, and the root-to-root order rule of the 10-24 exit applies (a move > 2 % on a per-watt figure ⇒ the clinical figures are quoted at degree 2 only); a gate red ⇒ a known-issues entry on the clinical degree-2 column, no adjudication; killed or OOM ⇒ the last flushed phase line is the lower bound, the modules split across two Saturdays. No band moves.

**Command:** written above by the 11-07 weekly (two calls, `timeout -k 60 14400` each, inside the 28 800 s window). No `XL_OUTPUTS`: the JSONs, the grid and the fields are gitignored.

### 63. `xl` (Friday 2026-11-13) — `GEO-36` step 10 / I32′: entry 55 re-run with the re-posed tip gate — the clinical implant-heating figure, degree 1, `-n 8` (weekly 2026-11-07, §10 chain I, entry 63; entry 55's family — the same mesh, drive, node, width and module, only the (S2-tip) comparand changed by I32b)

**Status:** **READ 2026-11-14 (weekly): RED, the branch ruled (PROJECT_PLAN §10 rulings 2026-11-14).** (a) The (S2-tip) red takes its pre-registered branch: the reading is recorded (R_Π 1.0920e-04 / 1.1039e-04 at Δt₁/10 against 1.1012e-03 / 1.0839e-03 at Δt₁, ratios 10.08 / 9.82 — the O(Δt) step term; the DG0 form gap is gone) and the tip read drops to printed on the clinical body; `S2_BAND` never moves. (b) The (E) red is not a chosen gate — I32b folded the Δt₁/10 BE step into the arm's asserted maximum (module `:354`, `:418–420`), a reading of the absolute-T floor by the 11-10 18:00 ruling's class (≈ 1.6e-08 at Δt₁/10 at the clinical source density) — so the asserted (E) returns to entry 55's definition (CN + the Δt₁ BE step) with the Δt₁/10 residual printed beside its floor estimate (§10 I32c, a slot step); `E_BAND` stays 1e-09 on the gates entry 55 passed. (c) The EM record and every θ record (equal to entry 55's at every printed digit) stand as records. **Re-run as entry 71 (I32″) once I32c lands;** known-issues 2026-11-06 stays open until 71 reads green; I76 is re-posed as I76′ on this log (the F2182 prints are untouched by the red gates); U153 and 70′ wait on 71 green; entry 67 rides I32c's gates (re-dated behind 71). Declined: a `src/` θ-carrying solve, widening either band, re-running this window unchanged. *Was:* **RUN 2026-11-13 — `20261113T080008Z_GEO-36-step10-I32-prime.log`: RED.** Status 1, `[capture] rc=1`, 995 s, `memory.peak` 47 214 108 672 B (43.97 GiB), `-n 8`, orphans 0 / 0; inside both readout brackets. One failing node, `test_geo36_i32_rod_heating`: `gates: failures ['gel twin (E) 9.862e-09', 'rod (E) 9.220e-09', 'rod (S2-tip) 1.104e-04']` (`:12532`). Every EM gate, (B), (E′), the global (S2) and (S0) are green; the printed (E) on the CN and Δt₁ steps are green (`:12514`, `:12529`) — the red (E) is the residual of I32b's extra BE step at Δt₁/10, folded into the arm's (E) and never printed alone (module `:354`, `:418–420`). The rod's R_Π at Δt₁/10 is 1.0920e-04 / 1.1039e-04 against 1.1012e-03 / 1.0839e-03 at Δt₁ (`:12521–12522`). **The pre-registered rule names neither a thermal (E) red nor this pair: the 11-14 weekly reads the branch** (mechanical read by the 2026-11-13 03:00 review from a `log-pathologist` reading; ledger row; known-issues 2026-11-06 "Entry 63's half"). No band moves. *Was:* **QUEUED 2026-11-08 for 2026-11-13** (18:00 review; `xl-queue.d/20261113-GEO-36-step10-I32-prime.env`, diffed identical to the command below). The command carries no sha to substitute, and I32b adds no knob: the clinical body passes `tip_projection=True` and `f2182=` itself (`test_geo36_rod_f_clinical_heating.py:589–593`); the default-off keywords bind only the other callers. The new print forms (the stub's F2182 / 1 g ball / tip-ball `<q>` / `<Pq>` / R_Pi / retired-form lines, `20261108T181211Z_GEO-36-I32b-L-W0b.log:105–114`) were cleared by the leak checker: `20261109T000435Z_LEAK-CHECK-I32-prime.log`, Status 0. *Was:* **READY 2026-11-08 (I32b landed, §9 item 890 landing commit)**. The 12:00 slot ran the stub proof on `main`, green at `-n 2` (`20261108T181211Z_GEO-36-I32b-L-W0b.log`), `-n 1` (`20261108T181245Z_GEO-36-I32b-L-W0.log`) and `-n 3` (`20261108T181312Z_GEO-36-I32b-L-W0n3.log`), and Wreg was green (`20261108T181335Z_GEO-36-I32b-L-Wreg.log`). The slot queued nothing. The next review substitutes the sha, checks the new printed forms against the leak guard (the tip-ball `<q>` / `<Pq>` / R_Pi / retired-form line and the F2182-shaped probe line) and writes the queue file. *Was:* **PENDING PREREQUISITE (written 2026-11-07 weekly):** §10 **I32b** — the (S2-tip) comparand re-posed in `tests/validation/test_geo36_rod_f_clinical_heating.py` and the thermal helper it shares with `WF-15` step 14 (the consistent-mass projection ⟨Πq⟩_B at the unmoved `S2_BAND` on one BE step at Δt₁/10; ⟨q⟩_B, the ratio ⟨Πq⟩_B / ⟨q⟩_B − 1 and the Δt₁ read printed; the F2182-shaped probe print, 1 mm beyond each tip at 900 s), stub-run green on I32a's stub with `OPS-83`'s checkerboard source (the new read passes where the DG0 form reads 1.6e-03). **The first review after I32b lands checks the module's new printed forms against the leak guard, marks READY and queues `xl-queue.d/20261113-GEO-36-step10-I32-prime.env`** (trailing-7-day count on 11-13, cutoff 11-07: 11-08, 09, 10, 11, 12 = five, six with it). If I32b has not landed by 11-12 18:00, the night goes to a priced-family variant under the daily's licence or unspent.

**Why:** entry 55 ran with every EM gate, (B), (E), (E′) and the global (S2) green and one red, `rod (S2-tip)`, which `OPS-83` (the stub) and `OPS-84` (I41's rod pair) showed to be the comparand's form — the DG0 ball mean against the scheme's consistent-mass projection, Δt-independent, no residual (known-issues 2026-11-06; §10 rulings 2026-11-07). The EM record stands; this window buys the thermal figure under a gate that compares the scheme with itself, so the first clinical-scale implant-heating figure can land as a record rather than a red.

**Price (measured, entry 55 — the same mesh, drive, node, width and module):** build 468.1 s, EM 521.1 s, thermal 12.9 + 7.0 s, window 1 015 s, `memory.peak` 44.12 GiB (`20261106T080008Z_GEO-36-step10-I32.log`); one extra BE step at Δt₁/10 per arm adds seconds. **Predicted ≈ 1 000–1 150 s, 42–50 GiB; readout brackets window [900, 1 800] s, `memory.peak` [40, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands as entry 55's — I36b's EM gates, (B), (E), (E′), the global (S2), (S0) — plus the re-posed (S2-tip) against ⟨Πq⟩_B at `S2_BAND` on both tip balls of the rod arm):** θ at the tip 1 g balls (ball means, stated) and at P37's six box probe points at t\* on both arms; the tip θ ratio arm M / arm R beside entry 55's 1.089 and I41's 1.10; ⟨Πq⟩_B / ⟨q⟩_B − 1 on each tip ball (the size of the form gap on the clinical mesh, beside `OPS-84`'s 6e-5 … 1e-2 on I41's pair); the F2182-shaped probe θ 1 mm beyond each tip at 900 s (printed); the rod's power share; per-phase time and RSS.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ **I32 lands as a record — the first clinical-scale implant-heating figure (no band on the θ ratio; no tissue, no compliance claim, no external check)** and known-issues 2026-11-06 retires on the clinical mesh; (S2-tip) red again with everything else green ⇒ the known-issues entry gains the reading and the tip read drops to printed — the comparand stays the question and `S2_BAND` never moves; an EM gate red ⇒ as entry 48's rule; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly.

**Command** (entry 55's command byte for byte except `XL_CHUNK` and the raw-log name; I32b changes the module, not the knobs — if I32b adds a knob, the READY review names it here and says so):

```
XL_CHUNK="GEO-36-step10-I32-prime"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo36 && R=/workspace/logs/geo36-i32-prime-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_HEAT=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the JSON and the fields are gitignored.

### 64. `xl` (Sunday 2026-11-15) — `GEO-38` C5f: three degree-2 nodes around the extrapolated root on F-clinical-3T at 128 MHz, `-n 16` (weekly 2026-11-09, §10 chain C, entry 64; entry 53's pattern on the 3 T twin — C12's module with the degree knob, entry 56's family for price)

**Status:** **RUN `20261115T080008Z_GEO-38-C5f.log`** (2026-11-15; marked by the 10:30 review). **READY and QUEUED 2026-11-09 (weekly) for Sunday 2026-11-15 02:00** — `docs/testing/xl-queue.d/20261115-GEO-38-C5f.env`, the command below, written with the Write tool. Prerequisites: C5e landed the 3 T degree-2 build path with (Deg) U₂ = 7 335 988 exact (`20261108T180244Z_GEO-38-C5e-W1-build.log`); C12's module copies `FEM_EM_GEO38_C12_DEGREE` into C4a's knobs under `1` (`test_geo38_f_clinical_3t_esr_root.py:24–27`, `:52–54`) and anchors (Deg-2) against C8's 3 T record; entry 57 read tonight (the degree-1 root this grid is centred from). Node forms (41.7500, 42.5000, 43.2500 and their `.12e` pF / F forms) cleared by the leak guard: `20261110T031116Z_LEAK-CHECK-C5F.log` (the first grid) and the re-check named in §10 tonight. Trailing-7-day count on 11-15 (cutoff 11-09): 11-09, 10, 11, 12, 13 = five, six with it. **The clerk diffs the file against this entry once.**

**Why:** entry 58's rule poses the 3 T `ANS` case (`ANS-12`, SPEC v1 tonight), and the 10-12 epitaph says a benchmark column is recorded at the production convention (degree 2) before it is adjudicated — so the 3 T twin needs its degree-2 root, as the 1.5 T twin needed entries 53 / 56 before `ANS-11` v2. The 11-07 weekly declined this window *as an answer to the asymmetric-class question* (that reading is C35's, entry 65); it is bought tonight as `ANS-12`'s serial prerequisite. C5e's 64×64 at degree 2 (an `xxl`) is **not** bought: the root is read on the ESR-root module's two feeds, as at 1.5 T.

**Grid (3 nodes, ≤ 4 decimals in pF): `41.7500, 42.5000, 43.2500`.** Centred on −3.5 % from the degree-1 root 44.0525 pF (entry 57), spanning −5.2 % … −1.8 %. The estimate: the 1.5 T root moved −4.05 % from degree 1 to degree 2 (193.5996 → 185.7629) and F-human's −2.30 % (118.2437 → 115.5209); the mechanism the records support is the lumped sheets' order-dependent stray reactance (the width law, known-issues 2026-10-11 retired 10-26), a near-field quantity that does not scale with ω, so the 3 T move is expected inside the two measured ones — ≈ 42.27 pF by the 1.5 T figure, ≈ 43.04 by F-human's. The grid leaves ≈ 0.5 pF below the lower estimate and ≈ 0.2 pF above the upper. A 0.75 pF spacing (1.7 %) is coarser than entry 53's 0.67 % in relative terms; the secant at degree 1 falls ≈ 5.6 % per 0.05 pF through the 3 T root (entry 57), so the chord over a 0.75 pF bracket will read high by more than entry 53's did — the refinement window reads the root at 0.25 pF (entry 56's pattern) and this window's job is the bracket and the degree-2 slope. The nodes do not collide with any entry-49 / 57 node (all above 43.9).

**Price (measured in the family; scaling stated):** entry 53 (three degree-2 nodes on the 1.5 T coil mesh at `-n 16`): build 368.2 s, Z50 pair 3 149.8 s, node pairs 2 970.4 / 3 258.1 / 3 582.0 s, window 13 351 s, `memory.peak` 480.47 GiB (`20261104T080008Z_GEO-38-C5b.log:13041`, `:13095–13096`). The 3 T coil mesh has 0.930× the degree-2 unknowns (7 335 988 / 7 885 674) and its build read 327.5 s (entry 57). Factorisation time scales at least linearly in unknowns: 330 + 0.93 × 12 960 ≈ 12 400 s linear, ≈ 11 500 s at N^1.5 ⇒ **predicted ≈ 11 500–12 500 s (80–87 % of the ceiling), 440–470 GiB; readout brackets window [9 000, 14 400] s, `memory.peak` [400, 512] GiB.** Timeout 14 400 s; 512 GiB limit (entry 53 used 93.8 %; 0.93 of that is 87 %). Entry 53's node pairs grew 2 970 → 3 582 s with the node; the same growth here is inside the prediction. **Runs past the 04:30 slot** (≈ 02:00–05:30), as entries 53 / 56 did with no symptom (`OPS-54`).

**Readout (record; C4a's gates asserted at their imported bands as in entries 49 / 57 — (A0) 1 057 895, census 64, (Deg) U₂ 7 335 988 exact, kinds and (P_RC) at each node, (P50) and NEG on the Z50 pair, (Z) at 128 MHz; (Z50) printed only):** Im Y_in at each node on all four series; the module's sign-change root on each series, **which the module never asserts, so a green window does not mean a root was found**; the degree-2 secant beside the degree-1 one (2.5–2.8e-2 S/pF, entry 57); the order move in pF and % beside the 1.5 T twin's −4.05 %; per-node times and `memory.peak`.

**Decision rule (pre-registered; the 11-16 weekly reads it):** every gate green and one − → + crossing on P33 (1, +1j) inside [41.75, 43.25] ⇒ C_in(F-clinical-3T, d2) is recorded at this grid's resolution, as a record, and **the 0.25 pF refinement (entry 56's two-node straddle) is posed by that weekly** as entry 68 for the first open night; every gate green and no crossing ⇒ the extrapolation is off by more than ≈ 0.4 pF — the measured degree-2 slope re-centres one more window, not a known-issues entry until the second grid misses; a gate red ⇒ a known-issues entry on the 3 T degree-2 fixture (entry 57's green at degree 1 is the comparison); killed at the timeout or OOM ⇒ the node JSONs written stand as records, the price is recorded, and the family moves to two nodes per window (entry 56's shape). No band moves. `ANS-12` SPEC v2 (the node) waits for the refinement; v1's record column is entry 58.

**Command** (entry 56's durable-capture shape on C12's module — entry 57's command with the degree knob added, the node list and the raw-log name changed; the clerk diffs it once):

```
XL_CHUNK="GEO-38-C5f"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c5f-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C12=1 FEM_EM_GEO38_C12_NODES_PF=41.7500,42.5000,43.2500 FEM_EM_GEO38_C12_DEGREE=2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_3t_esr_root.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the node JSONs are gitignored.

### 65. `xl` (Monday 2026-11-16) — `GEO-38` C35w: the degree-1 64×64 of F-clinical-1.5T on the r2 mesh, 64 MHz, `-n 16` (weekly 2026-11-09, §10 chain C, entry 65; entry 44's family — the same module, degree and width, the rung knob set; the h-reading of the asymmetric class)

**Status:** **RUN `20261116T080008Z_GEO-38-C35w.log`** (02:00 2026-11-16; Status 0, `[capture] rc=0`, 1 244 s; ledger row `xl-ledger.md:85`, filled by the 03:00 review; the decision rule is the 11-16 weekly's, with C51's class table, §9 item 1003). *Was:* **READY and QUEUED 2026-11-09 (weekly) for Monday 2026-11-16 02:00** — `docs/testing/xl-queue.d/20261116-GEO-38-C35w.env`, the command below, written with the Write tool. Prerequisite: C35 landed (§9 item 889) — `FEM_EM_GEO38_C3_RUNG=r2` on C3's module rebinds the kwargs and records to r2 and the save path to `logs/geo38/c3-cu-s64-d1-r2.npz` (`test_geo38_f_clinical_matrix.py:24–28`, `:64`, `:70`), build-only green at `-n 2` with 912 804 cells rel 0, U₁ 1 184 518 exact, census 64 (`20261108T193325Z_GEO-38-C35-W1-build.log`). No node, so no leak-guard form. Trailing-7-day count on 11-16 (cutoff 11-10): 11-10, 11, 12, 13, 15 = five, six with it. **The clerk diffs the file against this entry once.**

**Why:** known-issues 2026-10-27 (owner: the weekly) — the asymmetric class (other ring, 1 step) survived the order move ×1.44 on the 1.5 T twin (entry 54 / C18), so by that entry's second branch it is the fixture's or the mesh's; what separates them is an h-rung at degree 1 (11-07 ruling). The r2 rung (conductor_resolution 2.4e-3 m, `GEO-38` step 1a's second rung) is coarser than r1 (1.6e-3): a class whose spread falls with finer h is the mesh's; one that does not is the geometry's (the leg-box corner, the shield can).

**Price (measured in the family):** entry 44 on r1 (1 147 366 cells) at `-n 16`, degree 1: 1 233 s, 46.85 GiB — build ≈ 320 s, the factorising column ≈ 200 s, 63 held drives ≈ 9 s each. r2 is 0.80× the cells (912 804) and built in 245 s against r1's 323 s at `-n 2` (`20261015T171259Z_GEO-38-s1a-W2-r2-a2.log`) ⇒ **predicted ≈ 900–1 100 s, 35–45 GiB; readout brackets window [600, 2 400] s, `memory.peak` [25, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands as entry 44's — (A0) 912 804 rel 0 against C35's record, U₁ 1 184 518 exact, census 64, reciprocity at `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, `TH-14`'s copper identity, the 34 classes at the module's 5 % gate; the 0.5 % `ADJACENT_SPREAD_BAND` is the readout's band of record (entry 44's ruling), read, never a red):** the 34 C32 × mirror class spreads on r2 beside r1's (the asymmetric class 0.9262 %, the median, the max); the C13 statistic (other ring 1 step / same ring 1 step, r1 0.888); P_coil / P_in; the saved `c3-cu-s64-d1-r2.npz` (C51's input); per-phase times and `memory.peak`.

**Decision rule (pre-registered; the 11-16 weekly reads it, with C51's class table):** every gate green ⇒ the r2 64×64 is a record and C51 reads the class — **the asymmetric class's spread on r2 ≥ 2× r1's (it falls with finer h) ⇒ the class is the mesh's**, known-issues 2026-10-27 retires with the h-statement and no geometry step is posed; **its spread on r2 within 2× of r1's (no fall with h) ⇒ the fixture's**, the entry stays open and the weekly poses a geometry step (the leg-box corner first, the shield can second), one knob per window; a reciprocity or copper-identity red ⇒ a known-issues entry on the r2 rung's Leontovich wall; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command** (entry 44's durable-capture shape on the `xl` profile and service, C3's module with the rung knob; the clerk diffs it once):

```
XL_CHUNK="GEO-38-C35w"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-c35w-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C3=1 FEM_EM_GEO38_C3_RUNG=r2 timeout -k 60 14400 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_matrix.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS` beyond the saved `.npz` under the gitignored `logs/geo38/`.

### 66. `xxl` (Saturday 2026-11-21) — `GEO-38` C19-3T: the 3 T degree-2 column set at the 3 T degree-2 root — C6-3T (the box-loaded 3 T coil, `|B₁⁺|` / 10 g) and I36-3T (the rod in the box) at degree 2, 128 MHz, `-n 16` (weekly 2026-11-09, §10 chain C, entry 66; entry 62's twin on the 3 T fixture — the column `ANS-12` adjudicates on)

**Status:** **(iii) MET 2026-11-14 (weekly): entry 61 read, the first branch (the 3 T rod column at degree 1 is a record). (i) still open — entry 68 (`xl` 11-18 at the earliest) read by the 11-18 weekly, which writes the two-call command with the node, marks READY and queues `xxl-queue.d/20261121-GEO-38-C19-3T.env`; if 68 has not been read by the 11-20 18:00 review, entry 69's `20261128-…` file runs on 11-21 instead (lexical order) and this entry takes 11-28.** **(ii) MET 2026-11-12 (12:00 slot, §9 item 929 landing commit):** the knobs are `FEM_EM_GEO38_C21_DEGREE=2` (new, C21's module) and `FEM_EM_GEO36_S11_DRIVE_DEGREE=2` (new, C32's module); K-66 — C32's foreign guard refuses `FEM_EM_GEO38_C21` and `FEM_EM_GEO38_C21_DEGREE`, so the command is **two pytest calls, each with its own environment prefix** (entry 62's shape); W1-build `20261112T180748Z_GEO-38-C52-W1-build.log:11936` — (Deg-2) unknowns 7783028 = 2E + 2F = record, exact (E 1508256, F 2383258), the first degree-2 space on C20's mesh. (i) and (iii) still open. *Was:* **PENDING PREREQUISITES (written 2026-11-09 weekly):** (i) **the 3 T degree-2 root at 0.25 pF** — entry 64 (Sunday 11-15) read by the 11-16 weekly, then its refinement (entry 68, entry 56's two-node straddle; the first open `xl` night after 11-16); (ii) **§10 C52** — the degree-2 knobs on C21's module (`FEM_EM_GEO38_C21_DEGREE=2`, C34's pattern) and on C32's (`FEM_EM_GEO36_S11_DRIVE_DEGREE=2`), each stub-run with (Deg) U₂ = 2E + 2F exact at build on C20's mesh; (iii) **entry 61 read** (the 3 T rod column at degree 1, Thursday 11-12; the 11-14 weekly). **The first weekly after all three writes the command in entry 62's two-call shape with C21's / C32's modules, the degree knobs and the node, marks READY, and queues it into `xxl-queue.d/20261121-GEO-38-C19-3T.env`** (one `xxl` per trailing seven days: 11-21 is seven days after 11-14). If (i) is missing by 11-20 18:00, the Saturday goes to a priced-family variant under the daily's licence or unspent.

**Why:** `ANS-12`'s adjudication column, by the 10-12 epitaph — the 3 T twin of entry 62. The rod's degree-2 reading rides the same window because the factorisation is the price.

**Price (measured in the family; scaling stated):** entry 62's prediction (≈ 12 000–14 000 s, 500–540 GiB on the 1.5 T box mesh, 1 217 025 cells) scaled by the 3 T box mesh's 0.927× cells (1 128 430) and the 0.93× degree-2 unknowns ⇒ **re-priced 2026-11-14 (weekly) from entry 62's measured window — both calls 10 902 s (call 1 4 301.9 s with the box-wide grid, call 2 6 590.9 s), `memory.peak` 486.57 GiB at `-n 16` — scaled by the 3 T box mesh's 0.927× cells and 0.93× degree-2 unknowns ⇒ predicted ≈ 9 500–11 000 s, 440–480 GiB; readout brackets window [7 000, 18 000] s, `memory.peak` [380, 650] GiB.** *(Was: predicted ≈ 11 000–13 000 s, 465–500 GiB; brackets [8 000, 21 000] s, [400, 700] GiB.)* Timeout 28 800 s; 754 GiB limit; an OOM splits the two modules across two Saturdays.

**Readout and decision rule:** entry 62's, on the 3 T rows (entry 58's gates for C6-3T, entry 61's for I36-3T, (Deg) exact, (P_RC) at the node); every gate green ⇒ `ANS-12`'s `COMPARISON.md` degree-2 rows are filled and the order rule of the 10-24 exit applies; a gate red ⇒ a known-issues entry on the 3 T degree-2 column, no adjudication; killed or OOM ⇒ the modules split. No band moves.

**Command:** written by the READY weekly in entry 62's two-call shape (`R=/workspace/logs/geo38-c19-3t-xxl-raw.log`, `-n 16`). No `XL_OUTPUTS`.

### 67. `xl` (the first open night after entry 61 is read and I85 lands; Tuesday 2026-11-17 at the earliest) — `GEO-36` step 12 / I32-3T: the clinical implant-heating figure at 128 MHz — arm M's SAR heats the ASTM gel box around the rod inside F-clinical-3T, degree 1, `-n 8` (weekly 2026-11-09, §10 chain I, entry 67; entry 63's twin on the 3 T fixture — entry 61's mesh, drive and node, `MAT-15`'s thermal route and I32b's re-posed tip gate)

**Status:** **Re-dated 2026-11-14 (weekly): PENDING on (iv) §10 I32c landed and (v) entry 71 (I32″) read green mechanically (the 03:00 review's ledger row) — then READY: the first review after both writes the command in entry 71's shape with I85's module (`tests/validation/test_geo36_rod_f_clinical_3t_heating.py`, knob `FEM_EM_GEO36_S12_HEAT=1`) and the node 44.0525, `-n 8`, checks I85's printed forms against the leak guard, and queues it for Thursday 2026-11-19 (after entry 68 on 11-18; trailing-7-day count on 11-19, cutoff 11-13: 11-13, 15, 16, 17, 18 = five, six with it). (iii) entry 63 read RED and its branch ruled (PROJECT_PLAN §10 2026-11-14): this entry rides I32c's repaired asserted set — (E) on CN + Δt₁, the Δt₁/10 residual and the tip rows printed — not I32b's. Its decision rule is entry 71's, at 128 MHz.** *Was:* **PENDING on entry 63 read only — entry 63 ran RED on 11-13 (its (E) at Δt₁/10 and the rod's (S2-tip); EM green); this entry carries the same thermal route and Δt₁/10 gates, so whether it is bought, and in what form, is the 11-14 weekly's (03:00 11-13 review).** — **(ii) MET 2026-11-12 (03:00 review): entry 61 green on the mechanical read** (`20261112T080008Z_GEO-36-step11-I36-3T.log`, Status 0, every asserted gate green; ledger row); (iii) waits for the 11-14 weekly's read of entry 63. *Was:* **PENDING on entries 61 / 63 only** — (i) **§10 I85 MET 2026-11-10 (§9 item 932 landing commit; module `tests/validation/test_geo36_rod_f_clinical_3t_heating.py`, gated test `test_geo36_i32_3t_rod_heating`, knob `FEM_EM_GEO36_S12_HEAT=1` plus the node `FEM_EM_GEO36_S8_DRIVE_C_PF` (required); degree 1 only — the S8 / S11 knobs refused by its foreign guard; `20261110T105744Z_GEO-36-I85-W0b.log`, `20261110T105845Z_GEO-36-I85-W0n1.log`)**; (ii) **entry 61 green** (Thursday 11-12; read by the 11-14 weekly — the EM record the thermal arm rides on); (iii) **entry 63 read** (Friday 11-13; the 11-14 weekly — the 1.5 T figure's branch, so the 3 T twin is not bought on a gate that just went red again). **The first review after all three writes the command in entry 63's shape with I85's module and knob and the node, marks READY and queues it** for the first open `xl` night (11-17 at the earliest; the trailing-7-day count is the clerk's). If (i) is missing by 11-16 18:00, the night goes to a priced-family variant under the daily's licence or unspent.

*Was:* **PENDING PREREQUISITES (written 2026-11-09 weekly):** (i) **§10 I85** — the 3 T heating body, a slot step: I32b's module composed on C32's 3 T drive route (C20's rod-in-box mesh, 128 MHz, the node knob 44.0525, both arms), the re-posed (S2-tip) gate and the F2182-shaped probe print as I32b's, stub-run on I32a's 12³ stub with `OPS-83`'s checkerboard source; (ii) **entry 61 green** (Thursday 11-12; read by the 11-14 weekly — the EM record the thermal arm rides on); (iii) **entry 63 read** (Friday 11-13; the 11-14 weekly — the 1.5 T figure's branch, so the 3 T twin is not bought on a gate that just went red again). **The first review after all three writes the command in entry 63's shape with I85's module and knob and the node, marks READY and queues it** for the first open `xl` night (11-17 at the earliest; the trailing-7-day count is the clerk's). If (i) is missing by 11-16 18:00, the night goes to a priced-family variant under the daily's licence or unspent.

**Why:** the operator's directive names the 3 T twin second; the implant-heating figure at 128 MHz beside the 64 MHz one (entry 63) is the frequency lever on the rod's heating per accepted watt — the reading C40 / I76 are built for. A record: no tissue, no compliance claim, no external check (W3).

**Price (measured in the family):** entry 55 (the 1.5 T twin, `-n 8`): build 468.1 s, EM 521.1 s, thermal 12.9 + 7.0 s, window 1 015 s, `memory.peak` 44.12 GiB; the 3 T rod mesh is 0.927× the cells and entry 58 measured its build at 448 s ⇒ **predicted ≈ 950–1 150 s, 40–50 GiB; readout brackets window [800, 1 800] s, `memory.peak` [35, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout and decision rule:** entry 63's, at 128 MHz (I36b's EM gates as entry 61's, (B), (E), (E′), the global (S2), (S0), the re-posed (S2-tip) against ⟨Πq⟩_B at `S2_BAND`); the tip 1 g θ ratio arm M / arm R beside entry 63's and I41's 1.10; the F2182-shaped probe θ at 900 s printed; the 3 T / 1.5 T ratio of the tip θ per accepted watt (C40's twin, a record). Every gate green ⇒ the first 128 MHz clinical implant-heating figure lands as a record; (S2-tip) red ⇒ as entry 63's rule; an EM red ⇒ as entry 61's. No band moves.

**Command:** written by the READY review in entry 63's shape (`R=/workspace/logs/geo36-i32-3t-xl-raw.log`, I85's module and knob, `-n 8`). No `XL_OUTPUTS`.

### 68. `xl` (the first open night after the 11-16 weekly reads entry 64; Wednesday 2026-11-18 at the earliest) — `GEO-38` C5g: the 3 T degree-2 root refined to 0.25 pF on F-clinical-3T at 128 MHz, `-n 16` (weekly 2026-11-11, §10 chain C, entry 68; entry 56's pattern on the 3 T twin — C12's module with the degree knob, two nodes straddling entry 64's crossing)

**Status:** **PENDING PREREQUISITE (written 2026-11-11 weekly):** (i) **entry 64 read** by the 11-16 weekly — one − → + crossing of Im Y_in on P33 (1, +1j) inside the 41.75 / 42.50 / 43.25 pF grid, or its chord if the crossing falls between the outer nodes; if 64 finds no crossing inside its grid, the 11-16 weekly re-centres the grid instead and this entry becomes that re-run. **The 11-16 weekly writes the two straddling nodes (entry 56's shape: the bracketing pair's chord ± 0.125 pF, four decimals), checks their printed `.12e` forms against the leak guard, writes the command in entry 64's shape, marks READY and queues it** into `xl-queue.d/20261118-GEO-38-C5g.env` (or the first open night after 67 if 67 is READY first; the trailing-7-day count is the clerk's). Its reading writes `ANS-12` SPEC v2 (the node) and makes entry 66's (i) MET.

**Why:** `ANS-12`'s node and entry 66's prerequisite (i): the 10-12 epitaph puts the adjudication column at degree 2, and a 0.75 pF grid reads the root to its chord only; the straddle reads it at 0.25 pF as entry 56 did for the 1.5 T twin.

**Price (measured in the family; scaling stated):** entry 64's prediction (≈ 11 500–12 500 s, 440–470 GiB for three nodes plus the Z50 pair) × 2/3 for two nodes plus the pair ⇒ **predicted ≈ 8 000–9 000 s, 440–470 GiB; readout brackets window [5 000, 14 400] s, `memory.peak` [350, 512] GiB** — re-priced by the READY weekly from entry 64's measured per-node time. Timeout 14 400 s; 512 GiB limit.

**Readout (record):** Im Y_in on all four series at both nodes; the crossing by the bracketing pair on P33 (1, +1j); the four-series spread; per-node time and RSS.

**Decision rule (pre-registered):** one − → + crossing on P33 (1, +1j) inside the pair ⇒ **C_in(F-clinical-3T, d2) at 0.25 pF is the node**, written to four decimals as `ANS-12` SPEC v2 (R_p ‖ C and `Z_th` recomputed at 128 MHz) and entry 66's (i) MET; no crossing inside the pair ⇒ the chord of entry 64 is the node and the miss is a record; a gate red ⇒ the step returns to the weekly. No band moves.

**Command:** written by the READY weekly in entry 64's shape (`R=/workspace/logs/geo38-c5g-raw.log`, C12's module, `FEM_EM_GEO38_C12_DEGREE=2`, the two nodes, `-n 16`). No `XL_OUTPUTS`.

### 69. `xxl` (Saturday 2026-11-28) — `GEO-38` P33-d2: the clinical box heated by the coil at degree 2 — P61's module at the 1.5 T degree-2 root 185.7629 pF, both senses, `-n 16` (weekly 2026-11-11, §10 chain P / C, entry 69; entry 59's family at the order of entry 62 — the column `ANS-13`'s 1.5 T half adjudicates on)

**Status:** **READY and QUEUED 2026-11-14 (weekly) for Saturday 2026-11-28 02:00** — `docs/testing/xxl-queue.d/20261128-GEO-38-P33-d2.env`, the command below, written with the Write tool (the clerk diffs it once). (i) P87 MET (§9 item 969, below); (ii) entry 62 read tonight (its Status above; PROJECT_PLAN §10). Node forms (185.7629, its `.12e` pF / F forms, the knob string) cleared by the leak guard: `20261115T031801Z_LEAK-CHECK-P33-D2.log`, Status 0 — the same node entry 62 cleared (`20261108T030938Z_LEAK-CHECK-C19.log`). One `xxl` per trailing seven days: 11-28 is seven days after 11-21. **Lexical-order note:** entry 66's `20261121-…` file, once the 11-18 weekly writes it, runs first on 11-21; if 66 is not READY by the 11-20 18:00 review the launcher takes this file on 11-21 instead and 66 runs 11-28 — nothing is lost either way, and the clerk says which happened. *Was:* **PENDING PREREQUISITES (written 2026-11-11 weekly):** (i) **§10 P87** — the degree-2 knob on P61's module (`FEM_EM_GEO38_C6_DEGREE=2` passed through to C6a's patch set) with the volume-mean θ(t\*) print, stub-run with (Deg) U₂ = 2E + 2F exact at build and (B) on the degree-2 source — after entry 63 runs (the module hold on `test_geo36_rod_f_clinical_heating.py`, which P61 imports); (ii) **entry 62 read** by the 11-14 weekly — its measured window prices this one. **The first weekly after both writes the command in entry 59's durable-capture shape with the degree knob and the degree-2 node, prices it from entry 62's measured EM half plus the two thermal arms, marks READY and queues it into `xxl-queue.d/20261128-GEO-38-P33-d2.env`** (one `xxl` per trailing seven days: 11-28 is seven days after 11-21). If (i) is missing by 11-27 18:00, the Saturday goes to a priced-family variant under the daily's licence or unspent. **(i) MET 2026-11-13 (04:30 slot, §9 item 969):**
- **Knob:** `FEM_EM_GEO38_C6_DEGREE=2` routed through P61's `_run_heat_box(…, degree=)`.
- **(K-69):** the shape `FEM_EM_GEO38_C6_HEAT=1 FEM_EM_GEO38_C6_HEAT_SENSE=both FEM_EM_GEO38_C6_DEGREE=2 FEM_EM_GEO36_S8_DRIVE_C_PF=<node>` is accepted alone.
- **Stub proof:** (Deg) 8327394 exact (`20261113T105114Z_GEO-38-P87-W0a.log:331`), with (B) and (V12) at `B_BAND`.
- **Wording:** §10 calls this prerequisite "(ii)"; this entry calls it (i).

**Why:** `ANS-13`'s adjudication column at 1.5 T (the 10-12 epitaph: a benchmark column is recorded at the production convention before it is adjudicated); the degree-1 → degree-2 move of θ at the seven probes is the order sensitivity `ANS-9` probe D taught (8 % there).

**Price (measured in the family; scaling stated):** entry 62's prediction for C6 at degree 2 on this mesh (≈ 12 000–14 000 s, 500–540 GiB for the grid call; the box-wide grid is not needed here, so the readout half is entry 52's ≈ 125 s) plus two thermal arms of ≈ 8 s each ⇒ **re-priced 2026-11-14 (weekly) from entry 62's measured call 1 — C6 at degree 2 on this mesh at `-n 16`: build 480.2 s, MUMPS analysis + factorisation 166.6 + 2 691.8 s, degree-2 readouts 794.9 s (the box-wide grid's 626.7 s is not run here) — plus two thermal arms (≈ 20 s each at degree 1 on entries 55 / 59; the degree-2 source projection adds seconds) ⇒ predicted ≈ 4 000–5 000 s, 460–500 GiB (call 1 ran under the two-call window's 486.57 GiB peak); readout brackets window [2 500, 9 000] s, `memory.peak` [380, 600] GiB.** *(Was: predicted ≈ 9 000–12 000 s, 480–540 GiB; brackets [6 000, 21 000] s, [400, 700] GiB.)* Timeout 28 800 s; 754 GiB limit.

**Readout (record; asserted at their imported bands — the degree-2 EM gates of entry 62's C6 call; (B), (E), (E′), the global (S2); (Deg)):** `ANS-13` rows 1–6 at degree 2 (θ at c\* and P0–P5 per accepted watt under both senses, the factor P_box / P_acc, the calorimetric slope, the surface heat at t\*, the volume-mean θ(t\*)); the degree-1 → degree-2 move per row beside entry 59; per-phase time and RSS.

**Decision rule (pre-registered; the first weekly after the window reads it):** every assertion green ⇒ the degree-2 column lands and `ANS-13`'s 1.5 T rows are adjudicated on it when the operator's column arrives; a move > 10 % on any probe between orders ⇒ a record and the clinical thermal figures are quoted at degree 2 only (the 10-24 exit's order rule); a thermal identity red with the EM gates green ⇒ a known-issues entry on the thermal route at degree 2, the EM record stands; an EM gate red ⇒ as entry 62's rule; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. No band moves.

**Command** (entry 59's durable-capture shape with the degree knob and the degree-2 node — the (K-69) shape; `--profile xxl`, `-n 16`, one call under the window's 28 800 s ceiling; the clerk diffs the queue file against it once):

```
XL_CHUNK="GEO-38-P33-d2"
XL_COMMAND="docker compose --profile xxl exec -T fem-em-solver-xxl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo38 && R=/workspace/logs/geo38-p33-d2-xxl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO38_C6_HEAT=1 FEM_EM_GEO38_C6_HEAT_SENSE=both FEM_EM_GEO38_C6_DEGREE=2 FEM_EM_GEO36_S8_DRIVE_C_PF=185.7629 timeout -k 60 28800 mpiexec -n 16 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo38_f_clinical_box_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the JSON and the fields are gitignored.

### 70. `xxl` (Saturday 2026-12-05 at the earliest) — `GEO-38` P33-3T-d2: the 3 T clinical box heated by the coil at degree 2 — C31's module at the 3 T degree-2 root, both senses, `-n 16` (weekly 2026-11-11, §10 chain P / C, entry 70; entry 60's family at the order of entry 66 — the column `ANS-13`'s 3 T half adjudicates on)

**Status:** **PENDING PREREQUISITES (written 2026-11-11 weekly):** (i) **§10 P87** (C31's module honours the degree knob, with the volume-mean print); (ii) **the 3 T degree-2 root** — entry 68 read; (iii) **entry 66 read** (its measured window prices this one). **The first weekly after all three writes the command in entry 60's shape with the degree knob and the node, marks READY and queues it into `xxl-queue.d/20261205-GEO-38-P33-3T-d2.env`** (one `xxl` per trailing seven days: 12-05 is seven days after 11-28). If (ii) is missing by 12-04 18:00, the Saturday goes to a priced-family variant under the daily's licence or unspent. **(i) MET 2026-11-13 (04:30 slot, §9 item 969):**
- **Knob:** `FEM_EM_GEO38_C21_DEGREE=2` (C52's) routed through C31's `_run_heat_box_3t(…, degree=)`.
- **(K-70):** the shape `FEM_EM_GEO38_C31=1 FEM_EM_GEO38_C21_DEGREE=2 FEM_EM_GEO38_C6_HEAT_SENSE=both FEM_EM_GEO36_S8_DRIVE_C_PF=<node>` is accepted alone.
- **Stub proof:** (Deg) 7783028 exact (`20261113T105156Z_GEO-38-P87-W0b.log:316`), with (B) at `B_BAND`.

**Why:** `ANS-13`'s adjudication column at 3 T; its node is `ANS-12` v2's.

**Price (measured in the family; scaling stated):** entry 66's prediction (≈ 11 000–13 000 s, 465–500 GiB) for the EM half plus two thermal arms ⇒ **predicted ≈ 8 500–11 500 s, 450–510 GiB; readout brackets window [6 000, 21 000] s, `memory.peak` [400, 700] GiB** — re-priced by the READY weekly from entry 66's measured window. Timeout 28 800 s; 754 GiB limit.

**Readout and decision rule:** entry 69's, on the 3 T fixture, beside entry 60; the 3 T / 1.5 T lever at degree 2 beside C39's 4.0425 is the reading. No band moves.

**Command:** written by the READY weekly in entry 60's shape with the degree knob (`R=/workspace/logs/geo38-p33-3t-d2-xxl-raw.log`, `-n 16`). No `XL_OUTPUTS`.

### 71. `xl` (Tuesday 2026-11-17 at the earliest) — `GEO-36` step 10 / I32″: entry 63 re-run with I32c's asserted set — the clinical implant-heating figure, degree 1, `-n 8` (weekly 2026-11-14, §10 chain I, entry 71; entry 63's family — the same mesh, drive, node, width and module; only the asserted set changes: (E) on CN + the Δt₁ BE step as entry 55's, the Δt₁/10 residual and the tip rows printed)

**Status:** **PENDING PREREQUISITE (written 2026-11-14 weekly):** §10 **I32c** — the gate repair in `tests/validation/test_geo36_rod_f_clinical_heating.py`: `out["E"]` the maximum over CN and the Δt₁ BE step (entry 55's assertion), the Δt₁/10 residual printed as `(E10)` beside its floor estimate ε·T∞·ρc·V_sub / (Δt·∫q); a Δt₁/100 one-step BE row added to the tip rows with the two- and three-point intercepts printed; (S2-tip) printed on the clinical body and still asserted at `S2_BAND` on the stub; the T-from-0 separating measurement printed on the stub; I85's stub literals re-recorded under the (1\*) licence; stub-proved at `-n 2 / 1 / 3` with I32b's literals reproduced except the lines the brief lists. **The first review after I32c lands checks the new printed forms (`(E10)`, the Δt₁/100 row, the intercepts, the T-from-0 line) against the leak guard, marks READY and queues `xl-queue.d/20261117-GEO-36-step10-I32-dprime.env`** (trailing-7-day count on 11-17, cutoff 11-11: 11-11, 12, 13, 15, 16 = five, six with it). If I32c has not landed by 11-16 18:00, the night goes to a priced-family variant under the daily's licence or unspent, and this entry takes the first open night after I32c (before entry 67, which waits on it).

**Why:** entry 63 bought every EM gate and every θ record green and went red on two one-step readings — one the pre-registered (S2-tip) branch, one a gate I32b folded into the asserted maximum below the route's measured floor (PROJECT_PLAN §10 rulings and epitaphs, 2026-11-14). This window buys the same figure under the asserted set entry 55 passed plus the printed readings, so the first clinical-scale implant-heating figure can land as a record rather than a red, and known-issues 2026-11-06 can retire on the clinical mesh.

**Price (measured, entry 63 — the same mesh, drive, node, width and module):** build 451.1 s, EM 507.7 s, thermal 17.6 / 13.2 s, window 990.0 s, `memory.peak` 43.97 GiB (`20261113T080008Z_GEO-36-step10-I32-prime.log:12531`, `:12720`); one extra BE step at Δt₁/100 per arm adds seconds. **Predicted ≈ 1 000–1 150 s, 42–50 GiB; readout brackets window [900, 1 800] s, `memory.peak` [40, 70] GiB.** Timeout 14 400 s; 512 GiB limit.

**Readout (record; asserted at their imported bands — I36b's EM gates, (B), (E) on CN and the Δt₁ BE step, (E′), the global (S2), (S0); printed, no band — `(E10)` beside its floor estimate, the tip rows R_Π at Δt₁ / Δt₁/10 / Δt₁/100 with the two- and three-point intercepts, ⟨Πq⟩/⟨q⟩ − 1, the F2182-shaped probe θ at 900 s):** θ at the tip 1 g balls and at P37's six probe points at t\* on both arms; the tip θ ratio arm M / arm R beside entry 63's 1.0890 / 1.0898 and I41's 1.10; the rod's power share; per-phase time and RSS.

**Decision rule (pre-registered; a branch for every asserted gate class):** every asserted gate green ⇒ **I32 lands as a record — the first clinical-scale implant-heating figure (no band on the θ ratio; no tissue, no compliance claim, no external check)**, known-issues 2026-11-06 retires on the clinical mesh, and entry 67 is READY on the mechanical read; (E) on CN / Δt₁, (E′) or the global (S2) red ⇒ the step returns to the weekly, no second window (entry 63 printed them at 3.0e-12 / 4.3e-10, 3.3e-14 / 9.0e-14 and 1.25e-05, so this branch is not predicted); (S0) or (B) red ⇒ a known-issues entry on the source mapping, the EM record stands; an EM gate red ⇒ as entry 48's rule; killed or OOM ⇒ the last flushed phase line is the lower bound, the step returns to the weekly. The readings (`(E10)`, the tip rows) can turn nothing red. No band moves.

**Command** (entry 63's command byte for byte except `XL_CHUNK` and the raw-log name; I32c changes the module, not the knobs — if I32c adds a knob, the READY review names it here and says so):

```
XL_CHUNK="GEO-36-step10-I32-dprime"
XL_COMMAND="docker compose --profile xl exec -T fem-em-solver-xl bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs /workspace/logs/geo36 && R=/workspace/logs/geo36-i32-dprime-xl-raw.log && { echo [orphans-before]; pgrep -c python3; true; } > \$R 2>&1; PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_HEAT=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996 timeout -k 60 14400 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_heating.py -v -s --tb=short >> \$R 2>&1; rc=\$?; { echo [orphans-after]; pgrep -c python3; true; } >> \$R 2>&1; echo \"[XL] memory.peak bytes:\" >> \$R; cat /sys/fs/cgroup/memory.peak >> \$R; echo \"[capture] rc=\$rc\" >> \$R; cat \$R; exit \$rc'"
```

No `XL_OUTPUTS`: the JSON and the fields are gitignored.

### Weekly rulings on RUN entries, 2026-11-14 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 61 (`GEO-36` I36-3T, `xl` 11-12):** read — the first branch:
  every asserted gate green ⇒ **I36-3T is a record — the rod in the
  clinical box at 128 MHz, degree 1: tip 1 g Q 4.748–4.884, η_rod(T1)
  3.407, rod share 7.43e-05; against 1.5 T Q ≈ 1.08–1.10×, η 0.935×, share
  1.99× (C40)**; the 3 T implant-heating window (67) stays posed, re-dated
  behind 71. The 34.23 GiB under the bracket is a price reading. No band
  moves.
- **Entry 62 (`GEO-38` C19, `xxl` 11-14):** read — the green branch:
  every asserted gate green on both calls ⇒ **the degree-2 clinical column
  exists; `ANS-11`'s `COMPARISON.md` degree-2 rows are filled (C37) and
  the order rule does not fire (per-watt moves −0.46 / −0.74 / −0.49 /
  −1.26 %, node change included)** — the degree-1 figures stand as records
  with the move beside them; the rod's degree-2 reading is C56's
  (predicted to fire the rule). Its measured window re-prices 66 and
  prices 69. No band moves.
- **Entry 63 (`GEO-36` I32′, `xl` 11-13):** read — **RED, the branch
  ruled:** (S2-tip) takes its pre-registered branch (reading recorded, the
  tip read drops to printed on the clinical body); the (E) red is the
  Δt₁/10 step I32b folded into the asserted maximum — a floor reading by
  the 11-10 18:00 ruling's class, never a chosen gate — so (E) returns to
  entry 55's definition with the residual printed (I32c); EM and θ records
  stand as entry 55's; **re-run as entry 71** once I32c lands;
  known-issues 2026-11-06 open until 71 reads green. No band moves.
  Lesson (epitaph): a decision rule names a branch for every asserted gate
  class.
- **READY and QUEUED tonight:** 69 (P33-d2, `xxl` 11-28; priced from
  entry 62's call 1 ≈ 4 000–5 000 s, 460–500 GiB; lexical order lets 66
  run first on 11-21 if READY, else 69 runs 11-21 and 66 on 11-28).
  **Pre-registered PENDING:** 71 (I32″, `xl` 11-17 at the earliest, on
  I32c). **Re-dated:** 67 (11-19 at the earliest, on I32c and 71's
  mechanical green). **Re-priced:** 66 (≈ 9 500–11 000 s, 440–480 GiB;
  (iii) MET). 68 / 70 unchanged; 70′ (the Ti rod) the daily's variant
  licence after 71 green. **Not spent:** nothing — every `xl` night 11-15
  … 11-20 has a named candidate (64, 65, 71, 68, 67, 70′), the last two
  conditional.

### Weekly rulings on RUN entries, 2026-11-11 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 59 (`GEO-38` P33, `xl` 11-10):** read — the first branch: every
  *executed* assertion green (the readout's hotspot-ball slope and "(S0)"
  were never in the module; the DG0 tip print is known-issues 2026-11-06's
  scope) ⇒ **P33 is a record — the first clinical-scale coil-driven heating
  figure, θ(t\*) at c\* 1.612e-03 K per accepted watt at 1.5 T** (wall
  probes ≈ 2.6–2.8× hotter; `cem43` per watt the 20 °C baseline, not a
  dose) **and the 1.5 T thermal `ANS` case is posed: `ANS-13` SPEC v1.**
  No band moves.
- **Entry 60 (`GEO-38` P33-3T, `xl` 11-11):** read — the first branch:
  every executed assertion green ⇒ **P33-3T is a record — 6.516e-03 K per
  accepted watt at c\* at 3 T; the 3 T / 1.5 T lever 4.0425 ccw / 4.5397
  cw (C39) against the SAR10g lever 4.10 and the pictures' 3.69** (the
  pictures' shape part was 1 by construction — epitaph; P89) **and the 3 T
  thermal case is posed beside the 1.5 T one** (`ANS-13`'s second column,
  node with `ANS-12` v2). The 32.75 GiB under the bracket is a price
  reading. No band moves.
- **The after-window licence (§10 ruling):** a step "serial on entry N
  green" is takeable from the first slot after the 03:00 review's
  mechanical read (Status 0, rc 0, orphans 0, every asserted gate green by
  the pathologist row) — C39 is the precedent; entry 67's (ii) is marked
  MET the same way.
- **Pre-registered PENDING tonight:** 68 (the 3 T degree-2 root refined,
  `xl` 11-18 at the earliest, on entry 64), 69 (P33-d2, `xxl` 11-28, on
  P87 and entry 62) and 70 (P33-3T-d2, `xxl` 12-05, on P87, 68 and 66).
  66 / 67 unchanged. **Not spent:** the `xl` nights 11-17 … 11-20 hold 67
  / 68 and the daily's variant licence (the Ti-rod variant after I87′ and
  entry 63).

### Weekly rulings on RUN entries, 2026-11-09 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 57 (`GEO-38` C16r, `xl` 11-08):** read — the first branch: every
  gate green, one − → + crossing on P33 (1, +1j) at 44.052487 pF by the
  bracketing pair, all four series inside the grid (spread 0.0339 pF);
  **C_in(F-clinical-3T, d1) = 44.0525 pF at 0.05 pF** — an internal
  statement at degree 1; entries 58 / 60 / 61 carry the node. No band
  moves.
- **Entry 58 (`GEO-38` C6-3T, `xl` 11-09):** read — the first branch:
  every assertion green ⇒ **the first 128 MHz clinical-scale `|B₁⁺|` /
  10 g record** at the degree-1 root (lattice mean 1.175e-08 T, SAR10g(c\*)
  1.748e-04 W/kg at x+105 mm dz+50 mm, 9.328e-03 W/kg per accepted watt);
  3 T / 1.5 T per accepted watt 1.411 / 4.101 / 3.694 (lattice `|B₁⁺|`,
  SAR10g(c\*), P_box / P_acc; C33). **The 3 T `ANS` case is posed:
  `ANS-12` SPEC v1** (node TBD — the 3 T degree-2 root). **Entry 60 READY
  and QUEUED** (11-11). No band moves.
- **Pre-registered and QUEUED tonight:** 64 (C5f, the 3 T degree-2 root
  grid, `xl` 11-15) and 65 (C35w, the r2 h-reading, `xl` 11-16), both at
  measured family prices. **Pre-registered PENDING:** 66 (C19-3T, `xxl`
  11-21) and 67 (I32-3T, `xl` 11-17 at the earliest).

### Weekly rulings on RUN entries, 2026-11-07 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 54 (`GEO-38` C5d, `xxl` 11-07):** read — the second branch:
  every gate green at degree 2, one class (other ring, 1 step) at
  0.6424 % against 0.9262 % at degree 1; the degree-2 1.5 T 64×64 is a
  record (`C_mode`(d2) 185.4048 pF), the fixture's symmetry stays the
  finding, no band moves. Known-issues 2026-10-27's ≥ 2× clause does not
  fire (×1.44; the median class ×1.36) — the asymmetry is not degree-1
  discretisation; the 3 T degree-2 window is **not** commissioned on it;
  the next reading is §10 C35 (the r2 h-rung at degree 1). The log stays
  leak-held (one MUMPS timing line).
- **Entry 56 (`GEO-38` C5b′, `xl` 11-05):** read — the first branch: one
  − → + crossing on P33 (1, +1j), 185.762904 pF by the bracketing pair;
  the node is written as 185.7629 pF (entry 51's convention); P41
  (1, −1j)'s crossing at 185.820370 pF is a record (spread 0.057 pF
  against 0.117 pF at degree 1). **Entry 62 READY** at it, its two-call
  command written; `ANS-11` SPEC v2 carries it.
- **Entry 55 (`GEO-36` I32, `xl` 11-06):** read — the pre-registered
  "slope identity red with the EM gates green" branch; the EM record
  stands, the two Readout corrections are ratified (a 1 g ball mean; no
  Im Y_in print), and the (S2-tip) gate is re-posed as option (b) — the
  consistent-mass projection ⟨Πq⟩_B at the unmoved `S2_BAND`, one BE step
  at Δt₁/10 — for I32 and U127 together. **Entry 63 pre-registered**
  (I32′, Friday 11-13, PENDING on I32b, at entry 55's measured price).

### Weekly rulings on RUN entries, 2026-11-04 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 52 (`GEO-38` C6, 11-03):** read — the first branch, with
  two corrections: the registered readout was the 18-point lattice mean
  (≈ 1.498e-08 T by review arithmetic on `:12289–12306`), the module
  printed the 400-centroid mean (2.090e-08 T) under that name — both
  are the record, each with its definition; c\* (SAR10g 1.379e-04 W/kg,
  2.274e-03 per accepted watt, at x −105 / dz −50 mm) is the lattice
  argmax, not a box maximum (C6b, a slot step, adds the box-wide
  reader). C6 a record at degree 1, an internal statement; C7 posed as
  `ANS-11` (SPEC v1 written; its node is entry 56's); **entry 58 priced
  and READY except its node; entry 59 priced and READY on P61b.**
- **Entry 53 (`GEO-38` C5b, 11-04):** read — the first branch; the
  degree-2 root bracketed at 1.25 pF (185.763346 pF by the chord,
  −4.05 % from degree 1); all four series cross — an expectation miss,
  recorded, not routed on; the price growth (7.7 % over) is the MUMPS
  factorisation's, the slot overlap a correlate (`OPS-54`'s row: no
  symptom). **Entry 56 ratified** on the 10:30 review's straddling grid;
  **entry 62 (C19, `xxl` 11-14) pre-registered** on 56 + C34.
- **Entries 60 / 61 pre-registered** (P33-3T, I36-3T) on their slot
  bodies C31 / C32 and on entries 58 / 57.

### Weekly rulings on RUN entries, 2026-11-02 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 48 (`GEO-36` I36, 11-02):** read — the first branch; every
  gate that ran green on both arms, box containment not asserted on
  this route (the pathologist's correction, carried onto entry 55);
  I36 a record (tip 1 g Q 4.34–4.46, η_rod(T1) 3.642, the rod's power
  share 3.7e-5, S perturbed at 1e-7); I32 posed at this price.
  **Entry 55 kept** as queued for Friday 11-06, its price ratified.
- **Entry 49 (`GEO-38` C16, 11-01):** read — the first branch;
  C_in(F-clinical-3T, d1) = 44.0615 pF by the bracketing pair at 0.28
  pF; the sign of (root − `C_mode`) differs from the 1.5 T twin's and is
  C26's to read (a reader step, not a window). **Entry 57 kept** for
  Sunday 11-08; its node is entry 58's.
- **Entries 52 / 53 / 54 / 56:** unchanged, as the 10-31 rulings left
  them; the 11-04 weekly reads 52 and poses C7 and prices 58 / 59.
- **Pre-registered tonight, PENDING:** 58 (C6-3T, `xl` Monday 11-09 at
  the earliest, on C21 + 57 + 52 read) and 59 (P33, `xl` 11-10 at the
  earliest, on 52 green + P61). Nothing unspent 11-03 … 11-08.

### Weekly rulings on RUN entries, 2026-10-31 (full text PROJECT_PLAN §10; no AED figure opened tonight; the 10-28 weekly was lost to the host outage)

- **Entry 47 (`GEO-38` C14, 10-27):** read — branch 2 as annotated 10-26;
  the 3 T 64×64 is reciprocal, passive and copper-identity green at
  degree 1 (record); the one asymmetric class is known-issues
  2026-10-27, whose retirement condition is written (entry 54's 1.5 T
  order reading, C18; then C5e).
- **Entry 51 (`GEO-38` C4r, 10-28):** read — branch 1; C_in(F-clinical-1.5T,
  d1) = 193.5996 pF at 0.25 pF is the record; the bracketing pair's
  interpolation is the reading on every root grid; the 193.5996 node in
  entries 48 / 50 / 52 ratified; the held log is the operator's.
- **Entry 50 (`GEO-38` C5a, `xxl` 10-31):** read — first branch; the
  clinical degree-2 column priced (472.42 GiB, node pair 3 056.9 s);
  tier: root grids `xl` at ≤ 3 nodes (OOM ⇒ `xxl`, no `xl` re-run), the
  column set and root-node column sets `xxl`. **Entries 53 / 54 kept.**
- **Entry 52 (C6):** sense, lattice and coverage ratified; QUEUED
  Tuesday 11-03. **Entries 49 / 48:** unchanged, Sunday / Monday.
- **Commissioned:** 56 (C5b′, `xl` Thursday 11-05, PENDING on 53) and 55
  (I32, `xl` Friday 11-06, PENDING on 48 + I32a). Nothing unspent 11-01 …
  11-07; the 10-29 / 10-30 nights were lost to the host.

### Weekly rulings on RUN entries, 2026-10-26 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 45 (`GEO-36` step 8 / I34, 10-25):** read — first branch; the
  rod-in-box column priced (773 s, 41.76 GiB, exponent 1.643); I36
  commissioned, its body I36b, its night Friday 10-30 (entry 48
  re-dated); the Readout's "containment asserted" corrected (not on the
  column path).
- **Entry 46 (`GEO-38` C4, 10-26):** read — first branch; C_in(F-clinical-1.5T,
  d1) exists at 193.5178 pF (interpolated, 0.77 % below `C_mode`); C5 / C6
  posed at it; entry 51 kept to refine it; entry 50's node is entry 51's
  crossing (fallback 193.5178).
- **Entry 47 (C14):** QUEUED, unchanged; its second branch now reads C13's
  statistic (known-issues 2026-10-24 retired as a record). **Entry 49:**
  PENDING on 47 + C15, unchanged. **Entry 50:** re-conditioned (above).
  **Entry 51:** kept. **Entry 52 (C6): PENDING on C6a + I36b, `xl` Sunday
  11-01.** Nothing unspent 10-27 … 11-01.

### Weekly rulings on RUN entries, 2026-10-24 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 41 (`WF-10` step 4 / H28′, 10-22):** read — second branch
  (SAR10g −6.30 % root to root), every gate green ⇒ the human-scale tuned
  figures at degree 2 only, at C_in(d2); **Phase 6 exited** on the stated
  scope; H31's carry ratified.
- **Entry 43 (`IMP-3` step 10 / I30, 10-23):** read — first branch
  (3.539570 in [3.48, 3.56]) ⇒ the four-mesh record; known-issues
  2026-10-14 retired by its clause; the interval read as weak (epitaph);
  no finer rung.
- **Entry 44 (`GEO-38` C3, 10-24):** read — every gate green; the
  readout's 0.5 % band is the band of record ⇒ branch 2, the "other ring,
  one step" class at 0.926 % is the fixture's symmetry finding
  (known-issues 2026-10-24); C4 proceeds; C13 reads the saved matrix.
- **Entries 45 / 46:** kept (10-25 queued; 10-26 PENDING on the 10-25
  03:00 review's grid from C10). **Entry 47 (C14 = C11 full): READY — `xl`
  Tuesday 10-27, `-n 16`, command and price written tonight;** the 10-25
  03:00 review queues it. **Entry 48 (I36): PENDING on entry 45;** the
  10-26 03:00 review prices and queues it (`xl` 10-28). **Entry 49 (C16,
  the 3 T root): PENDING on 47 + C15;** the 10-28 03:00 review writes the
  grid (`xl` 10-29). **Entry 50 (C5a, the degree-2 cost probe): PENDING
  on entry 46's reading;** the 10-27 03:00 review sets the node (`xxl`
  10-31). **Not spent:** `xl` 10-30 (the 10-26 weekly fills it).
  **Backlog floor:** `xl` five ahead on the literal count (45, 46, 47,
  48, 49) — met; `xxl` one (50, PENDING) — met on the count. Declined: a
  finer tip rung; a degree-2 or h-refined 64×64 before C13 and C5a read.

### Weekly rulings on RUN entries, 2026-10-21 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 37 (`TH-20` step 2b / H13b, 10-20):** read — (A) 1.739e-11
  against 1e-9, inside both brackets ⇒ `TH-20` ✅ (the 10-20 03:00
  close-out ratified); the printed degree-2 controls get no window.
- **Entry 40 (`GEO-38` step 1b / C1 (b), 10-21):** read — every identity
  green, t_b 317.0 s, t_col 199.0 s, t_held 2.03 s, 40.37 GiB at `-n 8`
  ⇒ `GEO-38` ✅ as a probe; C3 priced from it; C5's degree-2 price
  re-read at the measured exponent 1.433 (`xxl`-only, after C4).
- **Entries 41 / 43:** kept as queued (10-22, 10-23). **Entry 44 (C3):
  READY — `xxl` Saturday 10-24, `-n 16`, command and price written
  tonight;** the 10-22 03:00 review queues it. **Entry 45 (I34): READY —
  `xl` Sunday 10-25, `-n 8`, command written tonight;** the same review
  queues it. **Entry 46 (C4):** PENDING C3 and C10; readout reconciled to
  C4a's module. **Backlog floor:** `xl` four ahead on the literal count
  (41, 43, 45, 46 PENDING) — met; `xxl` one (44) — met. Declined: a
  control window for entry 37's printed figures; re-running 37 or 40 for
  their price.

### Weekly rulings on RUN entries, 2026-10-19 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 36 (`WF-12` step 5 / H16c, 10-18):** read — gates green, the
  degree-2 copper column set exists; H25 read `C_mode`(d2, Cu) = 114.5221 pF
  from it, −0.865 % of the root (the circuit layer's prediction on copper
  at degree 2 inside 1 % of the in-model root); its share / leak
  predictions missed, printed, recorded not re-posed.
- **Entry 42 (`MAT-15` step 6 / P22 at the root, 10-19):** read — first
  branch: identities green, \|Im Y_in\| ≤ 1.91e-4 S, the ball-mean θ(t\*)
  per accepted watt −1.27 % root to root ⇒ P17's thermal figures quotable
  at degree 1 with the move beside them; `cem43` UNCOUNTABLE from the
  module → §10 P40 on the saved series; the price OUTSIDE-below is a
  record.
- **Entries 37 / 40 / 41 / 43:** kept as queued (10-20, 10-21, 10-22,
  10-23). **Entry 44:** PENDING entry 40 only; otherwise `xxl` 10-24 is not
  spent. **Added:** 45 (I34, `xl` Sunday 10-25, PENDING entry 40 and I35),
  46 (C4, `xl` Monday 10-26, PENDING C3 and C4a). **Backlog floor:** six
  `xl` entries ahead (10-20 … 10-26, two PENDING) — met; `xxl` one
  conditional entry. Declined: a Palace twin of entry 36; re-running 36 or
  42 for their price; a fourth entry-42 normaliser.

### Weekly rulings on RUN entries, 2026-10-17 (full text PROJECT_PLAN §10; no AED figure opened tonight)

- **Entry 34 (`ANS-6` step 5c / F31, 10-15):** read — the copper-loss
  partition converges (0.4760 → 0.4611 → 0.4576 on h = 0.015 / 0.0075 /
  0.005, ratio 0.23 on steps where first order gives 0.33); with T30's Q /
  `R_in` sequence (ratios 0.717 / 0.705 on steps where first order gives
  1.00) the h-statement is written and `ANS-7` is ✅ with the absolute Q
  INCONCLUSIVE (known-issues 2026-10-17; T48 asks the operator for AED's
  pass history). T30's "≤ ½" test is retired (§10 epitaph).
- **Entry 35 (`WF-10` step 4 / H20, 10-16):** read — second branch; the
  move is the detuning at the degree-1 node; the tuned figures are quoted
  at the degree-2 root only (entry 41).
- **Entry 39 (`WF-10` step 1e / H26, `xxl` 10-17):** read — first branch;
  C_in(d2) = 115.5209 pF; entry 41 READY and re-tiered to `xl` 10-22;
  `ANS-8` step 2's SPEC moves to this root (§10 H29) before the operator's
  session.
- **Entries 36 / 37 / 40:** kept as queued (10-18, 10-20, 10-21; 40 ratified).
  **Entry 38:** superseded by **42** (P22 at the degree-2 root; PENDING P22b;
  the 10-19 / 10-20 swap rule). **Entry 41:** READY, `xl` Thursday 10-22.
  **Added:** 43 (I30, `xl` Friday 10-23, READY), 44 (C3, `xxl` Saturday
  10-24, PENDING entry 40 and C3a — otherwise `xxl` 10-24 is not spent).
  **Backlog floor:** six `xl` entries ahead (10-18 … 10-23) — met; `xxl` one
  conditional entry. Declined: a Palace twin of entry 34, a fourth copper
  rung for Q, T43's h = 0.005 twin and a third `C_in` rung (standing).

### Weekly rulings on RUN entries, 2026-10-14 (full text PROJECT_PLAN §10; AED figures read for `ANS-9` only, private)

- **Entry 32 (`WF-12` step 4, 10-13):** read — the degree-2 copper F-human
  column is 514 s / 89.3 GiB at `-n 16` (drive 287 s), so by its own rule the
  full 32×32 is an `xl` case: entries 35–38 are commissioned on this price.
- **Entry 33 (`IMP-3` I15 r/8, 10-14):** read — Δ₁ +0.0739, Δ₂ −0.0258,
  opposite signs: the "otherwise" branch. The pointwise tip ratio is not
  resolved by the cap ball alone (known-issues 2026-10-14); §10 I17 opens; no
  finer cap-ball rung is bought (§10 epitaph).
- **Entries 34 / 35:** kept as queued (10-15, 10-16). **Entry 36:** `xl`,
  READY, Sunday 10-18, command written. **Entry 38:** READY, Monday 10-19 on
  its own night, command written (not folded into 35). **Entry 37:**
  re-worded (the mirror-built degree-2 column in the same window; step 1's
  module with the degree knob), READY, Tuesday 10-20, command written.
  **Added:** 39 (H26, the degree-2 root; PENDING H26a and entry 35's second
  branch; `xxl` 10-17 if READY by the 10-16 03:00 review, else `xl` 10-21).
  **Backlog floor:** five `xl` entries ahead (10-15, 16, 18, 19, 20) — met;
  `xxl` one conditional entry — **`xxl` 10-17 is not spent unless 39 is
  READY by Friday's 03:00 review**, and that is stated here rather than
  filled with an unpriced window. Declined: a Palace twin of T43 at
  h = 0.005 (no tuned node; not an `xl`-service job), `ANS-7` at 127.74 MHz
  on h = 0.005 (no AED column), the F-clinical degree-1 column before C1
  prices it.

### Weekly rulings on RUN entries, 2026-10-12 (full text PROJECT_PLAN §10; AED figures private)

- **Entry 30 (`ANS-4` 6c, 10-11):** read — F27's first branch: the self-class
  0.0075 → 0.005 change is 0.014 / 0.057 / 0.066 % at 10 / 64 / 128 MHz, so §2
  carries "h-converged to ≤ 0.1 % on three rungs"; the rung stays within 1 % of
  AED First Order on every class (the 10 MHz self class is the worst and does
  not move with h). No known-issues entry.
- **Entry 31 (`ANS-6` 5b, 10-12):** read — both clauses met; `ANS-6` row (a)
  AGREE within 0.33 % (×1, three frequencies) / 0.21 % (h = 0.0075, 64 MHz) of
  AED First Order on every class, the h-move toward it; `ANS-6` ✅. Entry 34
  (F31) stays queued as T30's copper-loss rung.
- **Entries 32–36:** kept as queued / pending. **Added:** 37 (`TH-20` 2b,
  PENDING entry 32; withdrawn if the quarter fits the heavy tier) and 38
  (`MAT-15` step 6 at degree 2, PENDING P22a and entry 32; may fold into 35).
  **Backlog floor:** three entries queued ahead (10-13 / 14 / 15) against the
  floor of four — unmet tonight, because every further candidate is serial on
  entry 32's price (10-13); the 10-14 weekly fills 10-18 / 10-19 from entries
  35–38 with measured prices. `xxl` 10-17: entry 36 or "not spent", decided
  10-14.

### Weekly rulings on RUN entries, 2026-10-10 (full text PROJECT_PLAN §10; AED figures private)

- **Entry 28 (`ANS-7` 4′, 10-08):** read — `f₀` AGREE with its h-statement
  (0.11 % at ×1, 0.03 % at h = 0.005 from AED First Order, converging first
  order in h; the rule's AGREE branch); `R_in` and Q move away (−2.7 % /
  −3.1 %), a stated level — T30 / entry 34 own the loss's h-sequence.
- **Entry 29 (`WF-7` 0e, 10-09):** not read — no `C_mode` exists for the saved
  matrix; H15c computes it (lane 1) and the 10-12 weekly reads rows 1 / 4.
- **Entries 30–33:** kept as queued. **Added:** 34 (READY, 10-15), 35
  (PENDING, 10-16), 36 (PENDING, `xxl` 10-17 or `xl` by entry 32's rule).

### Weekly rulings on RUN entries, 2026-10-07 (full text PROJECT_PLAN §10; no AED figure read tonight)

- **Entries 25 / 26 (`ANS-4` 4c / 4d, 10-06 / 10-07)** — decided: on the
  ragged family the degree-2 self-class conductor move falls with frequency,
  +1.92 / +0.48 / +0.24 % at 10 / 64 / 128 MHz, about 0.2× the first global
  h-step at Larmor, so **no conductor ladder is owed at 64 / 128 MHz**; the
  ragged family is closed (these were its last two windows, as the 10-05
  ruling said).
- **Kept:** 28 (10-08), 29 (10-09), 30 (10-11), 31 (10-12). **Added:** 32
  (H16 at `-n 16`; PENDING H16a; Tuesday 10-13) and 33 (I15 r/8 at `-n 16`;
  PENDING item 305's module; Wednesday 10-14). Both nights hold the
  six-per-trailing-7-days budget.
- **`xxl` 2026-10-10: not spent** — no candidate needs 754 GiB: the copper
  F-human degree-2 32×32 is priced by entry 32 and the σ = 800 degree-2 32×32
  by entry 29 first; the 10-10 weekly decides 10-17.

### Weekly rulings on RUN entries, 2026-10-05 (full text PROJECT_PLAN §10; AED figures private)

- **Entry 21 (`POST-7` 2, 10-04)** — decided: the ≈ 1 % degree-2 terminal
  residual at human scale is **not** an energy leak (exact identity
  1.207e-15 at 1e-6, NEG 431×). The weekly ruled a pre-registered cron window
  with logged assertions and a pathologist-read ledger row to be verification
  in §4's sense; `POST-7` closes on H7′ (the residuals registered as records).
- **Entry 27 (`ANS-4` 6b, 10-05)** — decided: T6 clause (i′) holds (self
  class moves 0.25 / 0.23 / 0.22 % against thresholds 2.32 / 0.94 %) and the
  rung is within 1 % of the First Order column on every class, as is the
  conductor-halved rung (items 212 / 225, clause (i″)). **T6 decided**:
  degree 2 on exact rectangles, κ off, is the convention; `PORT-21` ✅. At
  10 MHz the h-move is away from the column, so entry 30 is commissioned for
  the convergence statement, not the verdict.
- **Kept:** 25 (10-06), 26 (10-07), 28 (10-08), 29 (10-09) — 25 / 26 are the
  last ragged windows (a conductor move transfers). **Added:** 30 (READY,
  10-11), 31 (PENDING F24, 10-12). The daily's backlog-floor licence does not
  extend any ragged family beyond 25 / 26.
- **`xxl` 2026-10-10: not spent** — the copper F-human at degree 2 is priced
  by H16 (one column, an ordinary heavy window) first; if it lands by 10-07,
  the 10-07 weekly commissions the full 32×32 as `xl` entry 32.

### Weekly rulings on RUN entries, 2026-10-03 (full text PROJECT_PLAN §10; AED figures private)

- **Entry 16 (`ANS-6` 2d, 10-01)** and **entry 23 (`xxl`, `PORT-21` 3b,
  10-03)** — ragged-sheet records overtaken by the operator's 2026-10-02
  ruling (exact rectangles, κ off). Neither was read against the private
  column; neither decides anything. Entry 23's Status 1 is a commissioning
  error (its knob test needs a degree-1 h = 0.005 reference the entry never
  ran); no reference window is commissioned. Its in-run κ, 2.13e-2 at 10 and
  64 MHz, is recorded: κ shrinks with h on the ragged sheet.
- **Entry 20 (`ANS-4` 4b, 10-02)** — "comparable" is ruled **yes**: +1.92 %
  on the self class under conductor halving is 0.28–0.55× the global steps
  and 2.6× the gap the exact rectangles leave at 10 MHz. The conductor is a
  live h-axis; F18a carries it on exact rectangles. The `0.0075:2` /
  `0.005:2` ragged rungs are **not** commissioned.
- **Kept:** 21 (10-04), 25 (10-06), 26 (10-07). **Superseded:** 22.
  **24:** runs 10-08 unless 28 is READY by the 10-07 18:00 review.
  **Added, pending prerequisites:** 27, 28, 29. Monday 10-05 goes to 27 or
  29, whichever is READY at Sunday's 18:00 review.
- **`xxl` 2026-10-10: not spent** — nothing priced needs the tier now that
  human-scale degree 2 is measured at 103 GiB.

### Weekly rulings on RUN entries, 2026-09-30 (full text PROJECT_PLAN §10; AED figures private)

- **Entries 12, 13, 15 (3h, 3i, 3j)** — the congruent cut is not a term in
  the h-statement at 10, 64 or 128 MHz; each ladder reproduced its cut-off
  twin and none converges against its gap under the 09-26 criterion (last
  change on `S₁₁` 3.56 / 1.49 / 0.63 % against 5.7 / 2.5 / 1.1 %).
  **The global-resolution ladder family is closed**: no further window on
  it, and the daily's licence does not extend it.
- **Entry 14 (`ANS-6` 2c)** — `|S(Cu) − S(PEC)|` agrees with AED First
  Order to 0.5 % on every class at 64 MHz; the absolute classes sit 3.1 /
  1.8 / 1.5 % out, the same on both columns. Row (a) waits for entry 16.
- **Entry 16** kept (10-01); **20** kept (10-02); **21** kept (10-04).
  **17, 18 superseded; 19 superseded by 23.** Entries 22–24 added.
  **`xxl` 2026-10-03: spent on entry 23.**

### Weekly rulings on RUN entries, 2026-09-26 (full text PROJECT_PLAN §10; AED figures private)

- **Entry 5 (3d, 64 MHz, cut on)** — the cut is not a degree-2 systematic at
  64 MHz; two-rung moves within 0.012 pp of the cut-off rung's.
- **Entry 6 (`WF-7` 0c) / entry 9 (`WF-7` 0d, `xxl`)** — the F-human 32×32
  is heavy-tier at degree 1 and one ordinary window at degree 2 (103 GiB);
  the degree-2 terminal residual is on every column ⇒ entry 21.
- **Entries 7 and 11 (3e, 3g)** — both ladders move monotonically toward
  the order-matched AED column; **neither is converged against its gap**
  (last step ≈ 0.6× the gap). The "ratio ≳ 1.5" rule is **withdrawn** — the
  steps are unequal and the ladder never refines the port neighbourhood.
  3e's by-the-letter "converged" is not banked. T6 does not flip.
- **Entry 8 (3f)** — collapse branch ratified: the congruent cut is the rung
  to quote at 10 MHz.
- **Entry 10 (`ANS-2` step 4)** — branch (a) ratified; the pointwise
  statement is the distance to AED (within 6.3 %, 4.3 % drive-averaged), not
  our own 4.5 % scatter.
- **Entries 12, 13, 15 (3h, 3i, 3j)** — **kept**, and read under §10's
  2026-09-26 criterion (last successive change ≥ 3× below the quantity it
  is compared with), not under the ratio rule their text carries. Entry 14
  kept.
- Entries 16–21 added. **`xxl` 2026-10-03: not spent** — nothing priced
  needs the tier, and a fifth `ANS-4` rung extrapolates to ≈ 600 GiB on a
  754 GiB VM; not commissioned on an extrapolation.

### Weekly rulings on RUN entries, 2026-09-19 (full text PROJECT_PLAN §10; AED figures private)

- **Entry 1 (`WF-7` 0b, `xxl`)** — branch **(a)**: 5.50 % is in F-small's
  64 MHz class; degree 2 is the human-scale order at 106.1 GiB / 348 s per
  drive. The degree-2 terminal power residual (1.104e-02) is an unattributed
  record; chain H step H3 asserts the exact discrete identity.
- **Entry 2 (`ANS-4` step 3, 64 MHz)** — the pre-registered private rule
  **fired on the self class**; couplings inside it. The 64 MHz "AGREE by
  mechanism" is withdrawn for `S₁₁`, kept for the coupling classes.
  Known-issues 2026-09-19, `PORT-21`.
- **Entry 3 (`ANS-4` step 3b, 128 MHz, congruent cut)** — rule **did not
  fire**: classes within step 2d's spreads of step 2d's (≤ 0.083 % against
  ≤ 0.142 %). The cut is not a degree-2 systematic; the 128 MHz AGREE stands.
- **Entry 4 (`ANS-4` step 3c, 10 MHz)** — rule **fired**: the 10 MHz move is
  in the Larmor moves' class, so the degree-1 10 MHz AGREE is re-read as
  partly coincidental. The self-class spread at 89 % of its band is read
  with entry 8 when it runs.
- Entries 5–9 (daily licence) **kept**; entry 10 stays READY; entry 11 added.
