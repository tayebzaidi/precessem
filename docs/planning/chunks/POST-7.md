# POST-7 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §7 row text moved here verbatim — 2026-10-04 03:00 daily review (row at 2 KB)

> 🟡 **Step 1 DONE 2026-09-28** (§9 item 58): knobs `FEM_EM_WF7_DEGREE` / `FEM_EM_WF7_IDENTITY_COLUMNS` landed verbatim; degree 1, `-n 8`, heavy. W1 (knobs unset, 399.8 s window, 11.25 GiB summed RSS): (v) max rel 2.588e-15 on all 32 (band 1e-6); P17 reproduction `P_src,exact` rel 2.41e-11, κ rel 6.0e-8 (band 1e-6, asserted); NEG (one undriven sheet omitted) 4.793e-04 = **479× the band** (printed). W2 (`DEGREE=1 IDENTITY_COLUMNS=1`, 279.7 s window / 314 s pytest, 11.03 GiB): (v) 0.0, reproduction as W1, NEG 481× **asserted ≥ 10×**; terminal residual 4.187e-03 beside κ 4.590991e-03. Logs `20260928T093759Z_POST-7-step1-W1.log`, `20260928T094524Z_POST-7-step1-W2.log`. Step 2 = `xl` entry 21 (READY). **Step 1 (slot half, heavy):** a degree knob on `WF-7` step 1's gate module (unset = the landed degree-1 run, reproduced at the imported reproduction band), and `PORT-16`'s exact discrete identity asserted on the driven column at the imported `DISCRETE_IDENTITY_RTOL` with the terminal residual and the in-run κ printed beside it — at degree 1 this reproduces H3. **Step 2 (`xl`, `xl-pending.md` entry 21):** the same assert at degree 2, one factorisation, one driven column.

### Step 2 — the `xl` window, 2026-10-04 02:00 (read by the 03:00 daily review)

`20261004T070008Z_POST-7-step2.log`, cron `xl`, `-n 16`, command = `xl-pending.md` entry 21. A `log-pathologist` ruled it COUNTABLE; the review re-read `:10960–10961`, `:21689–21693`.
- 13 passed, 0 failed, 0 skipped on all 16 ranks; `[capture] rc=0` last; Status 0; 1097 s; `memory.peak` 104.39 GiB; orphans 0 / 0.
- Degree 2, 3 256 418 unknowns (`:10783`), factorisation 303.8 s; P17 drive 448.96 s; mesh 212.05 s; orientation-control rebuild 237.68 s.
- **(v) exact discrete identity on P17: max rel 1.207e-15, asserted at `DISCRETE_IDENTITY_RTOL` 1e-6** (`:10960`). NEG (P35's dissipation omitted) 4.311e-04 = 431× the band, asserted ≥ 10× (`:10961`; the printed ≥ 460× prediction missed). Terminal residual 1.104e-02 and pooled κ 1.215037e-02, printed.
- **Entry 21's decision rule: identity-green ⇒ degree 2 is not convicted;** the ≈ 1 % terminal residual is registered as a record at both orders (degree 1 4.187e-03 / κ 4.590991e-03, W2; degree 2 1.104e-02 / κ 1.215e-02). The κ-class reading is what the rule licenses; the log does not show residual = κ (κ is ≈ 10 % above it).
- Price: time OUTSIDE, above the 450–900 s bracket; memory INSIDE 100–115 GiB.
- One column only; (ii)–(iv) skipped by design under `IDENTITY_COLUMNS=1`. A self-consistency identity, not an accuracy or h statement.
- **Status:** the row's Done-when ("the identity asserted at degree 2 by an executed window, elapsed and `memory.peak` recorded") is met by this window. No chunk has been flipped to ✅ on a cron window before, and §4.2 says the agent executes the verification itself (`plan-navigator`, NOT FOUND for any rule or precedent). The review therefore holds 🟡 and names the ruling for the weekly; it does not set the precedent itself.

### §9 item 269 — 2026-10-06 (chain H7′, lane 1)

*Drafted for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 H7′ and the 2026-10-05 weekly's ruling "`POST-7` — §4.2 on a cron window, ruled" (PROJECT_PLAN.md §10, quoted, not changed): "`POST-7` closes when the identity-green branch's record registration lands (H7′, lane 1: the degree-1 and degree-2 terminal residuals as records in the `WF-7` module, no band, reproduced from `20260926T070009Z_WF-7-step0d.log` and `20261004T070008Z_POST-7-step2.log`); the ✅ is then the implementer's, audited as any other." No band, tolerance or validation target moves.

`POST-7` step 3: register the per-column terminal accounting residual at both element orders as records. No solve, no mesh.
- **Executor:** implementer. Tier **standard** (smoke-class content; the module's collection is the cost), `-n 2`, complex build. **Time:** ≈ 7 measured / ≈ 15 predicted min.
- **Depends on nothing.** The three logs are tracked, gzipped, under `docs/testing/logs/`.

**Correction to the licence's log citation (read before coding).** The degree-1 range is **not** in the step 0d log; step 0d is the degree-2 `xxl` run. The degree-1 range 4.061–4.624e-3 is step 0c's (`xl` 2026-09-21, `20260921T070008Z_WF-7-step0c.log`, `degree 1` at `:37`), as `xxl-ledger.md:21` itself says ("degree 1 at step 0c"). The records are read from the logs that printed them:

| Record | Value | Log : line |
| --- | --- | --- |
| degree 1, min over 32 columns | 4.061e-03 (drive 15) | `20260921T070008Z_WF-7-step0c.log:11187` |
| degree 1, max over 32 columns | 4.624e-03 (drive 20) | `20260921T070008Z_WF-7-step0c.log:11397` |
| degree 1, P17 | 4.187e-03 (drive 1) | `20260921T070008Z_WF-7-step0c.log:10599`; also `20260928T094524Z_POST-7-step1-W2.log` |
| degree 2, min over 32 columns | 1.081e-02 (drive 15) | `20260926T070009Z_WF-7-step0d.log:11187` |
| degree 2, max over 32 columns | 1.179e-02 (drive 20) | `20260926T070009Z_WF-7-step0d.log:11397` |
| degree 2, P17 (identity window) | 1.104e-02 | `20261004T070008Z_POST-7-step2.log:10955`, `:10960` |

Each line has the form `[WF-7 step0] phase 3 drive N accounting (printed, not asserted): … residual X.XXXe-0Y …` (step 0c / 0d), or `[WF-7 step1] drive 1/1 P17: … terminal residual 1.104e-02 …` (step 2).

**Code.** In `tests/validation/test_wf7_f_human_ring_matrix.py`, additive only:
- Module constants `TERMINAL_RESIDUAL_RECORD_DEG1 = (4.061e-03, 4.624e-03)`, `TERMINAL_RESIDUAL_RECORD_DEG2 = (1.081e-02, 1.179e-02)`, `TERMINAL_RESIDUAL_P17_RECORD = {1: 4.187e-03, 2: 1.104e-02}`, each with its log:line in a comment. Records, never bands: nothing in the module compares a fresh residual against them.
- One test, `test_post7_terminal_residual_records`, **not** `complex_only`, no mesh, no MPI collective beyond pytest's: rank 0 opens the three `.log.gz` files with `gzip.open(…, "rt")`, collects every `phase 3 drive N accounting (printed` line's `residual` token (exactly 32 per step 0 log — assert the count), and the step-2 `terminal residual` token. Broadcast the parsed values.

**Anchor (asserted, no band).** The registered constants equal the log prints **at print precision**: `f"{v:.3e}"` of each constant is string-equal to the min / max token of its log (and the P17 token), so no tolerance exists to be tighter than ½ ulp of a printed value. Plus two record identities, exact: each order's P17 value lies inside that order's 32-column range (4.061 ≤ 4.187 ≤ 4.624 e-3; 1.081 ≤ 1.104 ≤ 1.179 e-2), and the step-2 P17 token equals step 0d's drive-1 token (`:10599`, the same column at the same order — the ledger reports 0b/0d reprint it).

**Negative control (asserted; prior same-fixture measurement).** The degree-1 constants parsed against the degree-2 log (step 0d) miss: step 0d's min 1.081e-02 against 4.061e-03 (`:11187` of each log) — a 2.66× miss, string-unequal. Assert the swapped comparison fails for both min and max, and that the two ranges are disjoint (max deg-1 4.624e-03 < min deg-2 1.081e-02). A parser that read the wrong file or the wrong token would pass the anchor only by coincidence; this control shows it cannot.

**Window** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first). Price: `WF-7-H15a-W0`, the same module collected for one env-only contract test at `-n 1`, 100 s (`docs/testing/test-results.md:3057`). `-n 2` adds nothing material.
- **W1** (≈ 100–130 s), the closing window:

  `scripts/testing/run_and_log.sh POST-7-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py::test_post7_terminal_residual_records -v -s --tb=short'"`
- Node-ID selection, never `-k` inside the quoted command. **Do not** run the whole module: its default test is the `-n 8` heavy F-human solve.

**Traps already paid for.** Logs are `.log.gz` on disk, not `.log` (the licence's names are the harness names). `-k a or b` splits inside the quoted command. A node-ID run still imports every module the WF-7 module imports (H15a W0: 100 s for an env-only test), so the 30 s smoke ceiling does not hold — standard. No edit after W1 (rule (i)). The ledger's label traps on step 0d (`xxl-ledger.md:21`) concern prices, not the residual lines.

**Scope.** Covered: the per-column terminal accounting residual of the F-human 32×32 at 64 MHz, σ = 800 S/m volumetric conductor, ragged sheets as recorded, degree 1 (step 0c) and degree 2 (step 0d, step 2). **Not claimed:** that the residual equals κ (step 2 shows κ ≈ 10 % above it), any h-statement, any band on the residual, copper (`WF-12`), any AED reading.

**Negative result.** Any token mismatch ⇒ the constant is wrong or the line cited is wrong: fix the constant to the printed token (a record, rule (1*)-style re-read; never a band move), re-run W1. A count ≠ 32 ⇒ the parse is wrong, not the log; fix and re-run. Two failed windows ⇒ park on `attempt/POST-7-step3-*` and journal.

**Status it can move.** Done-when (§7 row `^| \`POST-7\` |`): "the identity asserted at degree 2 by an executed window, elapsed and `memory.peak` recorded. Green ⇒ … registered as a record at both orders". Step 2 met the first sentence (`:10960–10961`, 1097 s, 104.39 GiB); this step is the "registered as a record at both orders" clause — **it matches**. On green, `POST-7` 🟡 → ✅ by the implementer under the 10-05 ruling, audited as any other; §10 H7′ closes; the row's history moves here.

**Outcome — 2026-10-06 15:00 slot: GREEN, `POST-7` → ✅.** `20261006T200124Z_POST-7-step3-W1.log`, implementer, `-n 2`, Status 0, 71 s, 13 passed per rank (12 environment + `test_post7_terminal_residual_records`). `:79` degree 1 (step 0c): 32 columns, min 4.061e-03, max 4.624e-03, P17 4.187e-03; `:80` degree 2 (step 0d): 32 columns, min 1.081e-02, max 1.179e-02, P17 1.104e-02, step-2 P17 1.104e-02 (= step 0d drive 1); every constant string-equal at `.3e`. `:81` NEG (asserted): degree-1 constants against step 0d miss at min (2.66×) and max, ranges disjoint. One correction to the brief: the step-2 log is a plain `.log` on disk (step 0c / 0d are `.log.gz`); the helper opens `.gz` when present, else the plain file. No edit after W1. §10 H7′ is the weekly's to mark closed.

### §7 row status text replaced at closure — 2026-10-06 15:00 slot (verbatim)

> 🟡 **Step 2 green (`xl` 2026-10-04, `20261004T070008Z_POST-7-step2.log:10960–10961`): exact identity 1.207e-15 at degree 2 asserted at 1e-6, NEG 431×, 1097 s, 104.39 GiB — the Done-when is met; the 2026-10-05 weekly ruled that a pre-registered cron window with its logged assertions and a pathologist-read ledger row is verification executed in §4's sense — ✅ when step 3 (§10 H7′) registers the degree-1 and degree-2 terminal residuals as records in the `WF-7` module, no band.** Step 1 ✅ 2026-09-28 (degree 1, item 58).
