# ANS-12 — the ASTM-box-loaded F-clinical-3T birdcage under the two-port quadrature drive at its in-model Larmor root, 128 MHz

*History file for the §7 `ANS-12` row (`OPS-46` layout): the row holds the
state line and the Done-when; everything else lives here.*

## Commissioned — 2026-11-09 weekly planning review

Posed by `xl` entry 58's pre-registered decision rule ("every assertion
green ⇒ … the 3 T `ANS` case is posed by that weekly beside C7"): entry 58
(`20261109T080008Z_GEO-38-C6-3T.log`, 2026-11-09 02:00) ran with every gate
green at the 3 T degree-1 root 44.0525 pF (entry 57) and is the record column.
The adjudication column is the 3 T degree-2 set (`xl-pending.md` entry 66,
`xxl` Saturday 2026-11-21, PENDING on the 3 T degree-2 root — entry 64 then
its 0.25 pF refinement — and on C52 / entry 61).

**SPEC v1** (`examples/ansys_benchmarks/ans12_f_clinical_3t_box_loaded_128MHz/SPEC.md`):
the geometry, materials, 64 ESR capacitor sheets, the two feeds and the rows
are `ANS-11`'s with the 3 T fixture's dimensions and 128 MHz; the node is TBD
(SPEC v2, the weekly that reads the refinement). The row-3 tie clause is
`ANS-11` v2.1's.

**Slot steps (§10 chain C):** C48 — the `ans:12` record page (imports C21's
module, parses entry 58's tracked log through C33's parser; (X) / (D) / (R)
anchors as the `ans:11` page's); C49 — `COMPARISON.md` from C36's script
parametrised on the 3 T log, the degree-2 and AED columns blank, row 3's tie
set listed. Both lane 1, no solve.

**Operator (dashboard W11):** build the geometry now; enter the capacitor only
when SPEC v2 names it.

### §9 item 930 — 2026-11-10 (§10 C48, chain C — rule (4), the 11-09 weekly's step; lane 1 → chain step (no solve))

*Written by the 2026-11-10 03:00 daily review from the 11-09 weekly's C48 (`PROJECT_PLAN.md` §10 chain C), under §9 rule (4). It does not count against rule (6). Family `ANS-12` page, attempt 1 of 4. **Executor: implementer, in-slot.** This follows the `ans:11` precedent: item 833, the `ans:11` record page, and its two follow-ups C36 / C43 (items 871 / 880) were all run by the implementer (`ANS-11.md` outcomes). `example-runner` is for `EX-*` chunks. The 2026-11-05 03:00 ruling in `ANS-11.md` binds by analogy: a record reader, no solve, no mesh. The weekly's target binds: "entry 58's log tokens and the SPEC v1 table".*

**Two readings of the step text against the files (execution detail, target unchanged).**
1. **"Every SPEC-quoted token is found once" means each token is matched once, to its row's parsed field.** It cannot mean "occurs once in the log", which is false by construction for four tokens:
   - the centroid mean is printed at `:12235` and `:12237`;
   - the c\* SAR is printed at `:12253` and `:12256`;
   - the cw Y_in at P33 is printed at `:12233` and `:12234`;
   - S₁₂ = S₂₁ appears twice on `:12228` (all in `20261109T080008Z_GEO-38-C6-3T.log`).
   The raw multiplicity is printed.
2. **"(R) byte-regeneration of the SPEC's record table" is the lattice table** (`SPEC.md:178–182`, the header, the separator and three rows of 6 × 2 six-digit entries). Every byte of it follows from the log. The quantity table (`:163–173`) carries prose and a "Log line" column, and **that column is off on three rows**:
   - S is cited `:12227`; it is `:12228`;
   - Y_in ccw is cited `:12233`; it is `:12232`;
   - Y_in cw is cited `:12234`; it is `:12233`.
   The header's price citation `:12260` is the `wrote` line; the price line is `:12261`. So the quantity table goes through (X) token by token, and the page prints each parsed line beside the SPEC's. **The mismatches are a SPEC transcription question for the weekly, printed, never asserted, never edited.**

**Code: new files only.** No SPEC edit, no edit to C21's, C33's or C6a's module or to the `ans:11` page, and no `src/` edit.
- **`examples/ansys_benchmarks/ans12_f_clinical_3t_box_loaded_128MHz/12_f_clinical_3t_box_loaded_128MHz.py`.** Pure functions, then `main()`, wrapped so that a raise calls `MPI.COMM_WORLD.Abort(1)` (`imp:26`).
  - **Imports by name, never re-implemented.**
    - C21's module (`tests/validation/test_geo38_f_clinical_3t_box_loaded.py`): through `c21._modules()`, take C6a's `_lattice`, `BRIEF_LATTICE_N` and `C6A_LATTICE_DZ_M`, and `t3.FREQUENCY_3T_HZ` (printed).
    - C33's parser (`tests/validation/test_geo38_c33_field_strength_reader.py`): `_load_page`, `find_log`, `side`, `check_side`, `ratios`, `_half_ulp_rel`, `ENTRY58_GLOB`, `ENTRY58_NODE_PF` and `ENTRY52_NODE_PF`.
    - Through `_load_page()`, the `ans:11` page's `parse_c6_log`, `derived`, `_read_log`, `_spec_num`, `SOURCE_LOG` and `CW`.
    - Never import a `test_*` name. Never call `_c33_body`, `print_table` or `custody_json`.
  - **Reads.**
    - Entry 58 is `find_log(ENTRY58_GLOB)`, which must return `"tracked"`.
    - Entry 52 is `SOURCE_LOG`, used only for row `:173`'s ratios.
    - One local regex is new, for the line the `ans:11` parser does not cover: `(C21) mean |B1+|_ccw: centroid mean … lattice mean … over the 18-point lattice` (`:12235`).
    - One local SPEC parser is new, `parse_spec_3t(path)`. By regex, never by line number, it reads the quantity rows (`:165–173`), the lattice rows (`:180–182`), the tie-set sentence (`:184–186`) and the header's gate tokens (`:158–161`). The `ans:11` page's `parse_spec` does not match ANS-12's row titles.
  - **Prints.**
    - Every `[PASS]` / `[FAIL]` line.
    - The 18 rows with all digits and per P_acc: |B₁⁺|/√P_acc and SAR10g/P_acc, from `derived`. This is the SPEC's "the `ans:12` page carries all digits and the per-P_acc columns" (`:175–176`).
    - The tie set with its percentages.
    - The SPEC-line vs parsed-line table.
    - (N-sense): the cw / ccw P_acc ratio, printed only.
  - **Writes nothing** except (R-cmp)'s check, which only reads.
  - **Rank safety:** `_load_page()` on every rank before the rank-0 parse, because the page import is collective (`EX-193`, deviation 2). Run at `-n 1`.
- **Same-stem guide `.md`**, with `ans:11`'s headings:
  - the run line `./scripts/run_examples.sh -e ans:12 -n 1 -t 180`;
  - what is parsed, from which lines;
  - "a degree-1 record at the degree-1 root 44.0525 pF, not the node to enter; the adjudication column is entry 66";
  - the tie set;
  - links to `SPEC.md`, and to `COMPARISON.md` **only if item 931 has landed** (otherwise name it as C49's, unlinked).
  - It quotes only closing-log digits. Write `≈`, never a tilde.
- **Setup figure:** a `cmp`-identical copy of `examples/meshing/figures/meshing_16_f_clinical_3t_coil_astm_box_setup.png` (62 827 B, `mesh:16`'s render of the 3 T coil with the box), saved as `figures/ans12_f_clinical_3t_box_loaded_128MHz_setup.png`. Disclose it as `EX-193` did. The rod sub-volume is gel-mapped (SPEC `:73–80`).
- **`tests/unit/test_ans12_record_page.py`** (the `ans:11` test's loader shape, under its own `sys.modules` name) carries the gates and controls below.

**Gates (asserted in the script's `[PASS]` lines and in the unit test).**
- **(X) Extraction.**
  - The 18 lattice labels equal `_lattice()` in order.
  - The node string equals `ENTRY58_NODE_PF` ("44.0525").
  - Every SPEC token of rows `:165–172` equals its parsed field by string, after `_spec_num`: S ×3, Y_in ccw ×2 and cw ×2, the lattice mean (`:12235`), the centroid mean / |B₁⁻| / Sense / CV, c\* (SAR, cw ball, r, az, and x / z from the label's lattice point), P_box / P_src / P_acc, SAR per P_acc, and the centroid mean per √P_acc.
  - Row `:172`'s lattice mean per √P_acc (8.581609e-08) and row `:173`'s four ratios equal `side(p58)` / `ratios(side(p52), side(p58))` at `.6e`. C33 printed them at `20261109T103557Z_GEO-38-C33-W1-r2.log:195–198`.
  - The 36 six-digit entries (`:180–182`) equal `f"{float(tok):.6e}"`.
  - The header's gate tokens each occur in the log's C6a block:
    - (R1) 1.893e-13 at `:12228`; (P) 5.120e-13 / 1.810e-12 and NEG 1.2375e-02 / 1.2827e-02 at `:12229–12230`; (P_RC) at `:12231`; coverage 2.975e-14 / 3.220e-15 at `:12257`; build 448.3, columns 194.1 / 13.4 and readouts 117.1 at `:12261`;
    - 780 s is the footer's `Elapsed (s): 780` (`:12350`);
    - 33.48 GiB is `f"{35948146688 / 2**30:.2f}"` from `:12345`.
  - Every SPEC token is matched to exactly one row field, and none is unmatched.
- **(D) Derived, at half a unit in the last printed digit.**
  - C33's own `check_side(p58, ENTRY58_NODE_PF, "3 T (entry 58)")`, every check green:
    - (D-Pacc), `_p_acc_spec` against P_acc at ≤ 1e-8;
    - (D-box), (D-sar) and (D-cen), each at C33's Σ half-unit bound;
    - argmax = c\*, the c\* string equals its lattice line's, S₁₂ = S₂₁, and (Node).
    - C33 read these green on the same log (`…C33-W1-r2.log:186–193`): (D-Pacc) rel 5.755e-10; (D-box) 1.324e-08 ≤ 2.565e-07; (D-sar) 3.825e-08 ≤ 5.416e-08; (D-cen) 2.584e-07 ≤ 4.388e-07. They are therefore *asserted* (rule (e)).
  - **(D-lat):** the mean of the 18 parsed tokens equals the printed lattice mean at `:12235` within `_half_ulp_rel` (printed) + the max over the inputs. Review arithmetic (unlogged) gives equality to 10 digits, 1.174932617e-08. The `(C21)` line's centroid half equals `means.b1p_ccw` by string.
  - **(D-latP):** the lattice mean per √P_acc against SPEC's 8.581609e-08, at the same kind of bound.
- **(R) Regeneration.**
  - The page renders the lattice table from the parsed tokens. **It equals `SPEC.md:178–182` byte for byte**, including the Unicode minus in the row and column labels.
  - **The tie set.** The lattice balls with SAR10g within 1 % of the lattice maximum (the 11-09 tie clause, SPEC `:150`) are {x+105mm dz+50mm, x+105mm dz−50mm}. Each is listed with `f"{100 * (s / smax - 1):.2f}"`, and the nearest outsider x−105mm dz−50mm likewise. The strings −0.82 and −1.97 equal the SPEC's (`:184–186`; review arithmetic −0.8170 / −1.9700).
- **(R-cmp), conditional.** If item 931 has landed (ANS-12's `COMPARISON.md` exists and the `ans:11` page has its `"ans12"` case), the page re-renders it through that case and asserts byte identity. Otherwise it prints "absent — item 931 (C49)". This makes Done-when (1)'s closing window `-e ans:12` cover `COMPARISON.md` once both items have landed.
- `docrefs` `dead=0`.

**Controls (asserted).**
- **(N-digit-log):** in a `tmp_path` copy of entry 58's log, change the 4th significant digit of the c\* lattice line's SAR (`:12253`). (D-cstar), (X) and (R) go red.
- **(N-drop):** delete one lattice line. `parse_c6_log` raises (17 ≠ `BRIEF_LATTICE_N`).
- **(N-doc):** a one-digit change in a tmp SPEC copy's lattice table fails (R).
- **(N-tie):** in a tmp log copy, replace x+105mm dz−50mm's SAR token `1.734192670e-04` with `1.714192670e-04` (−1.96 %). The tie set has one member.
- **Printed only (rule (e)):** (N-sense), the cw / ccw P_acc ratio.

**Windows.** Bare, foreground, Bash `timeout: 660000`. `pgrep -c python3` in the container reads 0 first. `FEM_EM_ANS_SKIP_PRIVATE=1` on every pytest call.
1. Census-pre: `scripts/testing/run_and_log.sh ANS-12-page-census-pre "./scripts/run_examples.sh --list"`. The expected count is 188 (`EX-193`'s census-post, `20261110T043418Z_EX-193-census-post.log`), with ans 9. Read it; never assume it.
2. Docrefs-pre: `scripts/testing/run_and_log.sh ANS-12-page-docrefs-pre "timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py"`.
3. Code, guide and figure.
4. **W1** `ANS-12-page-W1`: `scripts/testing/run_and_log.sh ANS-12-page-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans12_record_page.py tests/unit/test_doc_reference_exit_codes.py -v -s --tb=short > /workspace/logs/ans12-page-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans12-page-w1-raw.log; cat /workspace/logs/ans12-page-w1-raw.log; exit \$rc'"`. This is C36's W1 form (`20261108T120606Z_ANS-11-C36-W1.log:12`; 39 passed, 23 s, `:977`, `:982`).
5. **Closing** `ANS-12-page`: `scripts/testing/run_and_log.sh ANS-12-page "./scripts/run_examples.sh -e ans:12 -n 1 -t 180"`. The `ans:11` page's closing window took 2 s (`20261108T152642Z_ANS-11-C43-page.log:126`). Nothing is edited after it, the guide included.
6. Census-post (189, ans 10), docrefs-post, and figures: `scripts/testing/run_and_log.sh ANS-12-page-figures "timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py"`.

**Tier.** Standard, `-n 1`, complex (the `ans:` group sources it; C21's import pulls in dolfinx), no solve, no mesh. **Class:** example (census +1; item 833's class). **Cost:** ≈ 10.3 measured / 50 predicted min. The windows sum to ≈ 1 min; the hour goes on the SPEC parser and the guide. If W1 is not green by minute 45, park.

**Traps already paid for.**
- Pytest output is interleaved on log lines: use `re.search`, never `startswith`. Retention gzips logs: read either suffix through `_read_log`, and cite the base name.
- **The `ans:11` page import is collective.** `_load_page()` on rank 0 alone hung `EX-193`'s first two W1s (`20261110T043240Z_EX-193-W1.log`).
- **C33's closing log holds two tables.** Its own self-test prints a table on a renamed copy (`EX-193`, deviation 1). This page computes the ratios through `ratios()` and never parses C33's log.
- **The leak hook** (`MIN_SIGNIFICANT_DIGITS = 8`): derived prints use ≤ 7 significant digits. Raw log tokens are quoted as the `ans:11` page quotes them. On a refusal, move the log unchanged to `logs/leak-held/`, report it, and never use `--no-verify`.
- Never write under `/workspace/logs/geo38` (C33's JSON; entry 66's outputs). Tmp copies go to `tmp_path`.
- Read no `aed_results/`, `COMPARISON_private.md` or `docs/private/`, and enter no AED value. The SPEC's status ticks are the weekly's.
- **`c21.K_NODE_3T_PF` is 44.0615** (entry 49's chord). Never print it as the node: the node is the log's (Z) token.
- Never delete a log, a red one included.

**Scope.** The degree-1 record page on entry 58. **Not claimed:**
- the degree-2 rows (entry 66);
- SPEC v2's node;
- an adjudication or any AED number;
- a box-wide maximum (the grid was off in entry 58);
- any tuning, compliance or homogeneity statement.

**Negative result.**
- **(X) red:** the parser or the log has moved, or a SPEC token is mistranscribed. Report the first mismatch, and edit neither file. A SPEC mismatch is the weekly's (SPEC v2). Mark BLOCKED.
- **(D) red with (X) green:** a printed figure does not follow from its printed inputs. That is a finding on C6a's writer. Park and mark BLOCKED; never retune.
- **(R) red:** the SPEC's lattice table is not the log's at six digits. Same disposition as (X).
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 C48 → DONE; census +1 (188 → 189, ans 9 → 10).
- **`ANS-12` ⬜ → 🟡 when items 930 and 931 have both landed.** That is Done-when (1): the page, (X) / (D) / (R), `COMPARISON.md` with the tie set, and the closing window. This brief reads the row's Done-when and matches it (rule (j)). Done-when (2) is the adjudication.

**Depends on:** nothing. Item 931 may run before or after it. If 931 has landed, (R-cmp) is live and the guide links `COMPARISON.md`.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 930 — `ANS-12` `ans:12` record page (C48) — per `docs/planning/chunks/ANS-12.md` "§9 item 930". New files only (page, guide, unit test, figure copy); no SPEC, C21, C33, C6a, `ans:11` or `src/` edit. No solve, no mesh. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (2026-11-10, 04:30 slot; implementer): DONE, every gate green.** Closing window `20261110T104306Z_ANS-12-page.log` (`-e ans:12 -n 1`, 2 s, exit 0): 61 `[PASS]`, 0 `[FAIL]` — (X) 46 lines (labels, node 44.0525, 28 quantity tokens, the 36 six-digit entries, 14 header gate tokens incl. 780 s / 33.48 GiB, the one-match coverage line "28 of 28, 14 of 14"); (D) C33's `check_side` 8 of 8 at C33's printed rel values, (D-lat) rel 1.419e-11 ≤ 7.679e-10 and the centroid string, (D-latP) rel 1.113e-08 ≤ 5.874e-08; (R) lattice table 662 B byte-identical to `SPEC.md:178–182`, tie set {x+105mm dz+50mm, x+105mm dz−50mm} at 0.00 / −0.82, outsider x−105mm dz−50mm −1.97. (R-cmp) printed "absent — item 931 (C49)". W1 `20261110T104239Z_ANS-12-page-W1-r2.log`: 33 passed, 23 s (the unit test's 9 nodes, controls (N-digit-log) / (N-drop) / (N-doc) / (N-tie) red as predicted). W1 r1 `20261110T104200Z_ANS-12-page-W1.log` was red 3 / 33 in `test_doc_reference_exit_codes.py` only: the copied setup PNG was not yet in the git index, so the pinned set (rule-admitted figure) exceeded `git ls-files`; staged, re-run green — no code change. Probe `20261110T104046Z_ANS-12-page-probe.log` (the script once, before W1). Census `…103758Z_ANS-12-page-census-pre.log` 188 (ans 9) → `…104314Z_ANS-12-page-census-post.log` 189 (ans 10); docrefs `dead=0` pre and post (exit 2, stale 121 unchanged, 188 → 189 guides); figures `…104318Z_ANS-12-page-figures.log` exit 0, `ans:12` ok. **Printed, never asserted (the weekly's, SPEC v2):** the SPEC "Log line" column is off on four citations — S `:12227` (parsed `:12228`), Y_in ccw `:12233` (`:12232`), Y_in cw `:12234` (`:12233`), header price `:12260` (`:12261`); the other five agree. (N-sense) cw / ccw P_acc ratio 1.0116, printed. **Deviations:** c\*'s r and az are compared as floats, not strings (the SPEC prints `0.105` / `0` beside the log's `0.105000` / `0.000`), and its x / z numerically against the label's lattice point, as the brief says; the item 931 helpers `TIE_REL` / `tie_set` are local to this page (931 adds its own to the `ans:11` page). `ANS-12` §7 row stays ⬜ until item 931 lands.

### §9 item 931 — 2026-11-10 (§10 C49, chain C — rule (4), the 11-09 weekly's step; lane 1 → chain step (no solve))

*Written by the 2026-11-10 03:00 daily review from the 11-09 weekly's C49 (`PROJECT_PLAN.md` §10 chain C), under §9 rule (4). It does not count against rule (6). Family C49, attempt 1 of 4. **Executor: implementer, in-slot** (C36 / item 871 precedent). The weekly's target binds: "the `ANS-11` rows regenerate byte-identical", "our degree-1 rows from entry 58, the degree-2 and AED columns blank until entry 66, row 3's tie set; `docrefs` 0 dead".*

**What "C36's comparison script" is.** It is not a separate file. Item 871 put it inside the `ans:11` page, `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/11_f_clinical_1p5t_box_loaded_64MHz.py`:
- `render_comparison` (`:667–737`), `build_metrics` (`:642–660`), `comparison_cells` (`:444`), `degree2_cells` (`:401`) and the `Q_*` labels (`:346–361`).
- "Parametrised, one branch" therefore edits that page. **C33's reader loads the page** through `_load_page()` (`test_geo38_c33_field_strength_reader.py:37–63`), and `mri:41` (`EX-193`) loads C33. Both are regression windows below.

**Code: the `ans:11` page, `tests/unit/test_ans11_record_page.py`, and the new `examples/ansys_benchmarks/ans12_f_clinical_3t_box_loaded_128MHz/COMPARISON.md`. Nothing else.**
1. **A `CASES` table.**
   - `"ans11"` holds today's literal strings: the title, the preamble, the column headers, Q_GRID's "not read in entry 52 (the C6b reader was off)", Q_CU's "… `:12282–12284`", and the provenance lines.
   - `"ans12"` holds the 3 T strings:
     - 128 MHz, entry 58, node 44.0525 pF "(entry 57's degree-1 root; not the node to enter, SPEC v2)";
     - the degree-2 header "ours, degree 2 (entry 66, pending)";
     - Q_GRID "not read in entry 58 (the C6b reader was off)" (no grid line in the log);
     - Q_CU "not printed by entry 58 (per-column partitions only, `:12229–12231`)";
     - the script name.
   - `render_comparison(m, c19=None, case="ans11")` and `build_metrics(parsed, d, case="ans11")`. **The defaults are byte-identical.**
2. **`SOURCE_LOG_3T_BASENAME = "20261109T080008Z_GEO-38-C6-3T.log"`**, in the shape of `SOURCE_LOG_BASENAME` (`:76`). **Do not import C33's module from this page:** C33 loads this page, so the import is a cycle.
3. **`tie_set(parsed, rel=TIE_REL)`**, with `TIE_REL = 0.01`, cited to the 11-09 weekly's tie clause (`ANS-11` v2.1; ANS-12 SPEC `:150`). It is the weekly's value, not a band set here. It returns the members (label, SAR string, `f"{100 * (s / smax - 1):.2f}"`) and the nearest outsider. In the `"ans12"` case it renders a "## Row 3 tie set" section; the `"ans11"` case renders none (the 1.5 T tie set is on the degree-2 column, C37).
4. **`main(argv=None)`** gains `--case ans12`. It parses entry 58, renders, writes ANS-12's `COMPARISON.md` only and prints. No ANS-12 `metrics.json` is written, so `LISTED_EXAMPLE_ARTIFACTS` is unchanged. The default path writes `ans:11`'s two files exactly as today. `run_examples.sh` passes no argument.
5. **No 3 T degree-2 reader.** The node is TBD (SPEC v2), so `c19=None` and the degree-2 column is blank.

**Gates (asserted in the unit test).**
- **(K) The `ANS-11` rows regenerate byte-identical.** The default render of entry 52 equals the committed `COMPARISON.md` (4 566 B), and `build_metrics` equals `metrics.json` (6 025 B). After the closing `-e ans:11` window, `git diff --stat` on both files is empty. Every existing `test_ans11_record_page.py` node is green; C43's W1 read 42 passed per rank (`20261108T152614Z_ANS-11-C43-W1-r2.log:1943`).
- **(X3T).** Every degree-1 cell of ANS-12's `COMPARISON.md`, read by `comparison_cells`, contains entry 58's parsed token string for its row: the 17 main rows (raw cells) and the 18 lattice rows (|B₁⁺|, SAR, and their log lines `:12238–12255`). Derived cells use the page's existing formats (≤ 7 significant digits).
- **(Tie).** The section lists exactly {x+105mm dz+50mm (0.00, the maximum), x+105mm dz−50mm (−0.82)}, plus the nearest outsider x−105mm dz−50mm (−1.97). The two percentages equal ANS-12 SPEC `:184–186`'s strings, parsed by regex.
- **(Blank).** The degree-2 and AED cells are `""` on every main row.
- **(R3T).** Re-rendering from the log equals the committed ANS-12 `COMPARISON.md` byte for byte.
- `docrefs` `dead=0` (`COMPARISON.md` links `SPEC.md`).

**Controls (asserted).**
- **(N-case):** entry 52 rendered with `case="ans12"` fails (X3T), and entry 58 rendered with `case="ans11"` differs from the committed ANS-12 file.
- **(N-tie):** in a `tmp_path` log copy with x+105mm dz−50mm's SAR token set to `1.714192670e-04`, the tie set has one member.
- **(N-doc):** a one-digit change to a tmp copy of ANS-12's `COMPARISON.md` differs from the render.
- C36's existing controls stay as they are.

**Windows.** Bare, foreground, Bash `timeout: 660000`. `pgrep -c python3` reads 0 first.
1. Code.
2. **Wgen** `ANS-12-C49-gen`: `scripts/testing/run_and_log.sh ANS-12-C49-gen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/11_f_clinical_1p5t_box_loaded_64MHz.py --case ans12'"`. It writes ANS-12's `COMPARISON.md`.
3. **W1** `ANS-12-C49-W1`: C36's W1 command (`20261108T120606Z_ANS-11-C36-W1.log:12`), labelled, with raw log `ans12-c49-w1-raw.log`. Add `tests/validation/test_geo38_c33_field_strength_reader.py` to the paths, knob unset, since it loads this page. Expect C43's 42 plus the new nodes plus C33's.
4. **Closing (K)** `ANS-12-C49-ans11-page`: `scripts/testing/run_and_log.sh ANS-12-C49-ans11-page "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"`. 71 `[PASS]` at C43 (`20261108T152642Z_ANS-11-C43-page.log`, 2 s, `:126`). Then run `git diff --stat` on `ans:11`'s `metrics.json` / `COMPARISON.md`: it must be empty.
5. **Wreg** `ANS-12-C49-mri41`: `scripts/testing/run_and_log.sh ANS-12-C49-mri41 "./scripts/run_examples.sh -e mri:41 -n 2 -t 120"`. This is `EX-193`'s closing form (`20261110T043404Z_EX-193.log`, 4 s, `:117`; 35 of 35 gates). It is a regression only and edits nothing.
6. **If item 930 has landed:** `scripts/testing/run_and_log.sh ANS-12-C49-ans12-page "./scripts/run_examples.sh -e ans:12 -n 1 -t 180"`. (R-cmp) is then live and green. This is the closing window of Done-when (1).
7. **Docrefs** `ANS-12-C49-docrefs`: item 871's form, `dead=0`.

**Tier.** Standard, `-n 1` (`mri:41` at its own `-n 2`), complex, no solve. **Class:** chain step (no solve). **Cost:** ≈ 5.0 measured / 35 predicted min. The windows sum to ≈ 1.5 min.

**Traps already paid for.**
- Item 871's list: `re.search`, never `startswith`; either log suffix; ≤ 7 significant digits on derived values; `FEM_EM_ANS_SKIP_PRIVATE=1`; a regenerated file is restored unless byte-identical.
- **No import cycle.** This page must not import C33's module.
- **No edit to `parse_c6_log`, `derived` or `parse_spec`:** C33 and `mri:41` read them.
- The page import is collective (`EX-193`): any test that loads it under `mpiexec` loads it on every rank.
- The leak hook: never `--no-verify`. On a refusal, move the log unchanged to `logs/leak-held/`.
- No AED file is opened. Never delete a log.
- **No module or page is edited after its last window** (rule (i)). If the page changes after the `-e ans:11` window, re-run it.

**Scope.** The comparison renderer's 3 T branch and ANS-12's `COMPARISON.md` with the degree-1 column. **Not claimed:**
- any degree-2 value or reader (entry 66; the node is SPEC v2's);
- the 1.5 T tie set (C37);
- an adjudication or any AED number;
- a SPEC edit.

**Negative result.**
- **(K) red:** the parametrisation leaked into the 1.5 T render. Fix the default path, keep the first log, and re-run once (rule (i)).
- **(X3T) red:** a parse or label error. Report it and stop, BLOCKED.
- **(Tie) red against the SPEC strings:** a SPEC transcription question for the weekly. Edit neither file; BLOCKED.
- **A C33 or `mri:41` regression red:** stop, BLOCKED, and leave the known-issues decision to the review.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 C49 → DONE.
- With item 930 landed, `ANS-12` ⬜ → 🟡 (Done-when (1)). Otherwise `ANS-12` is "one page (C48) away", as the weekly wrote.
- `ANS-11` stays 🟡.

**Depends on:** nothing. It may run before item 930; it edits none of 930's files.

**Outcome (2026-11-10, 04:30 slot; implementer): DONE, every gate green.** The `ans:11` page gains `CASES` (`"ans11"` today's strings, `"ans12"` the 3 T strings), `SOURCE_LOG_3T_BASENAME`, `COMPARISON_MD_3T`, `TIE_REL = 0.01`, `tie_set`, `render_comparison(m, c19=None, case="ans11")`, `build_metrics(parsed, d, case="ans11")` and `main(argv=None)` with `--case ans12` (writes ANS-12's `COMPARISON.md` only; no ANS-12 `metrics.json`); `parse_c6_log` / `derived` / `parse_spec` untouched, no C33 import. Wgen `20261110T104926Z_ANS-12-C49-gen.log` (2 s, exit 0) wrote 5 051 B and printed the tie set x+105mm dz+50mm (0.00), x+105mm dz−50mm (−0.82); outsider x−105mm dz−50mm (−1.97). W1 `20261110T104946Z_ANS-12-C49-W1.log` (`-n 2`, 23 s): **67 passed, 1 skipped per rank** (C43's 42 + 9 new C49 nodes + C33's file + `test_ans12_record_page.py`; the skip is C33's knob-gated `test_geo38_c33_reader`, knob unset as briefed). Gates: (K) default render 4 566 B and `build_metrics` 6 025 B equal the committed files; (X3T) 0 misses over 17 main rows + 18 lattice rows (tokens and `:12238–12255` line cells); (Tie) members / outsider and the SPEC `:184–186` percentage strings −0.82 / −1.97 by regex; (Blank) 17 rows, degree-2 and AED `""`; (R3T) 5 051 B byte-identical. Controls: (N-case) entry 52 as `ans12` gives 45 (X3T) misses, entry 58 as `ans11` differs; (N-tie) one member; (N-doc) differs. Closing (K) `20261110T105022Z_ANS-12-C49-ans11-page.log` 71 [PASS], 2 s, then `git diff --stat` on `ans:11`'s `metrics.json` / `COMPARISON.md` empty. `mri:41` `20261110T105032Z_ANS-12-C49-mri41.log` 35 of 35, 2 s. `ans:12` `20261110T105035Z_ANS-12-C49-ans12-page.log` 62 [PASS], (R-cmp) live 5 051 B vs 5 051 B, 2 s. Docrefs `20261110T105043Z_ANS-12-C49-docrefs.log` `dead=0`, stale 121 unchanged (exit 2). `LISTED_EXAMPLE_ARTIFACTS` unchanged: `.md` is not in `ARTIFACT_SUFFIXES`. **Deviations:** none to a gate; the tie section is a 4-column table (point, SAR string, %, role); 930's guide still names `COMPARISON.md` unlinked (930's file, not edited here). `ANS-12` ⬜ → 🟡 (Done-when (1)); §10 C49 left to the weekly.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 931 — `ANS-12` C49, `COMPARISON.md` — per `docs/planning/chunks/ANS-12.md` "§9 item 931". Edit only the `ans:11` page script and `tests/unit/test_ans11_record_page.py`, and write ANS-12's `COMPARISON.md` through the script's `--case ans12`; `ans:11`'s `metrics.json` / `COMPARISON.md` must stay byte-identical. No solve. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

### §9 item 989 — 2026-11-15 (§10 C58, chain C — rule (4), the 11-14 weekly's step; chain step (no solve))

*Enumerated by the 11-14 weekly at `PROJECT_PLAN.md:5211–5219`. Brief written by the 2026-11-15 03:00 daily review as execution detail. A record: no band is set or moved. Family C58, attempt 1 of 4. **Executor:** the implementer, in-slot (items 930 / 983 precedent). The weekly's target binds: "the entry-66 reader — C37's 3 T twin written ahead of its window: the `ans:12` page's degree-2 column parses entry 66's tracked log when the ledger names it, with a tree-independent "(Pend)" premise (the item 918 pattern) while it does not; stub-tested on entry 62's log shape — (K) the same parser reproduces C37's 53 degree-2 fields from `20261114T080009Z_GEO-38-C19.log` by string; `docrefs` 0 dead. Validation target: (K) exact; no solve."*

**What the step is.** C37 (item 983) gave the `ans:11` page a degree-2 column parsed from `xxl` entry 62's tracked log. C58 writes the 3 T twin's reader into the **`ans:12` page** before `xxl` entry 66 (C19-3T) runs. The reader finds entry 66's log through `docs/testing/xxl-ledger.md`. While the ledger names no such row, the page prints "(Pend)", and its premise test proves that state on empty tmp roots, not on the real tree. (K) runs the same reader on entry 62's tracked log and requires C37's 53 degree-2 field strings, exactly.

**Facts verified against the tree (this review, by `grep` / reading; unlogged, not a record).**
- **The `ans:12` page** is `examples/ansys_benchmarks/ans12_f_clinical_3t_box_loaded_128MHz/12_f_clinical_3t_box_loaded_128MHz.py` (503 lines).
  - It loads the `ans:11` page through `c33._load_page()` (`:71–73`, a collective import) and reads entry 58 through that page's `parse_c6_log` (`:221`).
  - `main()` (`:448–489`) prints the `[PASS]` / `[FAIL]` lines, and any `[FAIL]` returns 1.
  - Last window: `20261114T181222Z_ANS-11-C37-ans12-page.log`, 3 s, **62 `[PASS]`**, 0 `[FAIL]`, (R-cmp) live 5051 B == 5051 B (`:102`).
  - Its unit test `tests/unit/test_ans12_record_page.py` has 9 nodes. They are green on both ranks in `20261114T181121Z_ANS-11-C37-W1.log`: 76 passed per rank over `test_complex_mode.py` + the `ans11` / `ans12` / `ans13` page tests + `test_doc_reference_exit_codes.py`, 24 s, command at `:12`.
- **ANS-12's `COMPARISON.md` is written by the `ans:11` page**, not by the `ans:12` page: `main_ans12`, `11_f_clinical_1p5t_box_loaded_64MHz.py:1108–1125`, through `CASES["ans12"]`. The `ans:12` page only re-renders it for (R-cmp) (`:412–425`).
- **C37's parser** (the `ans:11` page):
  - `parse_c6_log` `:146–172` (any second match of a singleton line raises);
  - `_GRID_RE` / `parse_c19_grid` `:748–766` (exactly one grid line, else `ValueError`);
  - `degree2_block` `:811–825`, which pins the node to `C19_NODE_PF = "185.7629"` (`:368`);
  - `custody_c19` `:898–942`.
- **The 53 fields** are `custody_c19`'s pair names: `s2.s11/s12/s21/s22` (4); `y_in.{ccw,cw}.{p33,p41}`, each `"<re> <im>"` (4); `lattice.<label>.{b1p_ccw,sar10g}` for 18 labels (36); `c_star.label/sar10g/sar10g_cw` (3); `mean_b_plus.ccw` (1); `p_box/p_src/p_acc` (3); `grid.star_label/sar10g` (2). C37 measured `n_fields 53, mismatches 0` against the host JSON (`ANS-11.md` § item 983 Outcome).
- **The tracked comparand for (K)** is `ans:11`'s `metrics.json`, key `degree2_c19`, written by C37's closing window `20261114T181245Z_ANS-11-C37-page-r.log`. Read at this review, it carries `source_log` `20261114T080009Z_GEO-38-C19.log`, plus `parsed` (S2 `:12379`, Y_in ccw `:12383`, the 18 lattice lines from `:12388`, c\*, means, normalisers) and `grid` (`:12407`, `x+157.5mm y+0.0mm dz+0.0mm`), all as printed strings.
- **Entry 66** (`docs/testing/xl-pending.md:2547–2557`):
  - chunk `GEO-38-C19-3T`, so its log is `<stamp>_GEO-38-C19-3T.log`;
  - the queue file `xxl-queue.d/20261121-GEO-38-C19-3T.env` is **not written**: only `20261128-GEO-38-P33-d2.env` is queued;
  - the 11-18 weekly writes the command in entry 62's two-call shape, after entry 68.
- **`docs/testing/xxl-ledger.md`** has its header at `:18` and rows at `:20–27`, with **no `GEO-38-C19-3T` row**.
  - The 1.5 T row `:27` has chunk `GEO-38-C19`, a prefix of the 3 T name: a glob or regex `GEO-38-C19*` matches both.
  - Row `:26` (C5d) shows a leak-held log's form: the backticked base name, then "— **held by the private-leak guard, not in git**".
- **Entry 66's node is not written.** It is `ANS-12` SPEC v2's (the 11-18 weekly's, after entry 68), so no node constant can be asserted. The reader prints the node token only.
- **The grid on entry 66 depends on its command.** Entry 62's call 1 set `FEM_EM_GEO38_C6_GRID=1` (`20261114T080009Z_GEO-38-C19.log:12`), and the log has one result line (`:12407`; `:12162` is the census line, which `_GRID_RE` does not match). Neither foreign list names the grid knob (`test_geo38_f_clinical_box_loaded.py:63–73`, `test_geo38_f_clinical_3t_box_loaded.py:67–75`), and item 819's (K8-grid-3T) ran C6b on the 3 T box (`:581`). Whether entry 66 carries it is the 11-18 weekly's choice.
- **Item 918's pattern** (`GEO-38.md` § "§9 item 918"): a premise test monkeypatches the reader's roots to empty `tmp_path` directories (rank 0 creates them, then a bcast) and asserts "pending" by construction. It never asserts absence on the real tree.
- **Premise corrected (1): where the degree-2 column can land.** The `ans:11` page directly imports `xxl` entry 69's own module: `import tests.validation.test_geo38_f_clinical_box_heating as p61` in `main_ans13`, `11_f_clinical_1p5t_box_loaded_64MHz.py:1315`; `tests/unit/test_ans11_record_page.py:723` does the same. **So both are held until entry 69 runs (11-28).**
  - C58 therefore **cannot** fill ANS-12's `COMPARISON.md` degree-2 cells, because `CASES["ans12"]` lives in that held page.
  - "The `ans:12` page's degree-2 column" is the `ans:12` page's reader and its printed column. ANS-12's `COMPARISON.md` stays byte-identical, and (R-cmp) checks that.
- **Premise corrected (2): what (K) compares.** The weekly's "53 fields" is (J62)'s count, a JSON-against-log custody list. (K) compares the new reader's strings with C37's committed parse (`metrics.json`), which is tracked. It never reads the gitignored host JSON.
- **Module hold (the 03:00 11-15 digest, `PROJECT_PLAN.md:3985–3989`).** C58 edits the `ans:12` page script, `tests/unit/test_ans12_record_page.py` and, optionally, the page guide `12_f_clinical_3t_box_loaded_128MHz.md`.
  - None is a queued window's own module, and none is in entry 64's, 65's or 69's forward import closure (AST walk at this review: 69 / 55 / 116 modules).
  - None directly imports any of those windows' own modules (grep, both directions; the only other loaders of the page script are its own test and two path constants, `SPEC_3T` / `COMPARISON_MD_3T`).
  - The page imports C33's module (free) and C21's `test_geo38_f_clinical_3t_box_loaded.py`. C21 is in entry 69's closure and is itself a direct importer of entry 65's module (`:773`), so it is held; importing it is allowed, and C58 does not edit it.
  - The page also loads the held `ans:11` page; again allowed, and not edited.
  - **Result: takeable now.**

**Code.** Two files, plus the guide optionally (one sentence, no digits). Every other file stays byte-identical, including the `ans:11` page, its test, C33, C21, both SPECs, both `COMPARISON.md` files and `ans:11`'s `metrics.json`.
1. **Constants (page):**
   - `XXL_LEDGER = _REPO_ROOT / "docs" / "testing" / "xxl-ledger.md"`;
   - `TRACKED_LOG_DIR = _REPO_ROOT / "docs" / "testing" / "logs"`;
   - `ENTRY66_CHUNK = "GEO-38-C19-3T"`;
   - `ENTRY62_LOG_BASENAME = "20261114T080009Z_GEO-38-C19.log"`;
   - `D2_FIELDS_N = 53`, commented as C37's (J62) count, a record;
   - `_ROW66_RE = re.compile(r"^\| \d{4}-\d{2}-\d{2} \| GEO-38-C19-3T \| `(\d{8}T\d{6}Z_GEO-38-C19-3T\.log)`")`. The chunk column is matched exactly between pipes, so row `:27` cannot match.
2. **`entry66_log(ledger=XXL_LEDGER, log_dir=TRACKED_LOG_DIR)`.** It collects the base names matched across the ledger and prints `(rows, distinct)`.
   - 0 distinct ⇒ `None`.
   - More than 1 distinct ⇒ `ValueError` naming the count.
   - Exactly 1 ⇒ the path under `log_dir`, either suffix. If neither suffix exists (a leak-held or untracked log), it raises `ValueError` naming the file: the page reads tracked logs only, and a held log is the review's.
3. **`degree2_fields(path)`** — a pure function with the `ans:11` page's parsers **unedited**:
   - it runs `p = page.parse_c6_log(path)` and counts `page._GRID_RE` matches per line (`re.search`);
   - 0 matches ⇒ no grid fields (51 fields, "grid: not run" printed); 1 ⇒ `page.parse_c19_grid(path)` (53 fields); more than 1 ⇒ `ValueError`;
   - it returns `{"log", "node": (token, line), "fields": [(name, token)], "grid_lines": n}`, with names exactly `custody_c19`'s, in its order.
4. **`fields_from_metrics(path=page.METRICS_JSON)`** builds the same `[(name, token)]` from `metrics.json["degree2_c19"]["parsed"]` / `["grid"]`. It reads a tracked file and nothing else.
5. **`check_k(…)`** returns `[(name, ok, detail)]`. Do not reuse the `ans:11` page's `check_k` name in a way that shadows it: call this one `check_k_d2`. It holds:
   - **(K-n)** `len(fields) == D2_FIELDS_N` on entry 62's log;
   - **(K)** every name's token `==` the metrics string (exact string equality), with mismatches listed;
   - **(K-src)** `metrics["degree2_c19"]["source_log"] == ENTRY62_LOG_BASENAME == page.C19_LOG_BASENAME`, and the parsed node token `== page.C19_NODE_PF`.
6. **`main()`:**
   - appends `check_k_d2`'s lines to `results`, so the page's `failed` covers them;
   - prints `[ANS-12 page] degree 2: (Pend) entry 66 pending — xxl-ledger.md names no GEO-38-C19-3T row (rows 0, distinct 0)` while `entry66_log()` is `None`;
   - once the ledger names entry 66, it prints that log's node token and line, its `grid_lines` count and the 53 / 51 fields with their lines, and asserts only that the parse succeeds.
   - Every existing print and `[PASS]` line stays unchanged, so the count goes 62 → 65.
7. **Unit test nodes** (rank 0 writes every tmp file and broadcasts the path — C45's `_on_root`, `test_geo36_c40_rod_frequency_reader.py:67–79`, own copy):
   - `test_c58_k` — (K-n) / (K) / (K-src) on the real tracked files. This is tree-independent, because both files are tracked and entry 62 has run.
   - `test_c58_pend` — **(Pend)**: a tmp ledger holding only the real header lines and an empty tmp log dir ⇒ `None`. A tmp ledger holding only a verbatim copy of row `:27` (1.5 T) ⇒ `None`.
   - `test_c58_stub66` — **(Stub66)**: a tmp ledger with one row `| 2099-01-01 | GEO-38-C19-3T | \`20990101T000000Z_GEO-38-C19-3T.log\` | 16 | … |` and a tmp dir holding entry 62's log under that name (`page._read_log`, written uncompressed). The reader returns 53 fields equal by string to `test_c58_k`'s, with node token `185.7629`.
   - Controls: `test_c58_n_twin`, `_n_absent`, `_n_digit`, `_n_grid2` and `_nogrid`, below.

**Gates (asserted; no band set or moved).**
- (K-n), (K) exact, (K-src);
- (Pend) on empty roots; (Stub66);
- every existing `ans:12` gate unchanged: (X), (D), (R), and (R-cmp) live with ANS-12's `COMPARISON.md` 5051 B byte-identical;
- every existing node in the W1 set green, unedited;
- `docrefs` `dead=0`.

**Negative controls (asserted).**
- **(N-twin):** a tmp ledger with two distinct `GEO-38-C19-3T` base names raises, naming 2. The same row repeated is accepted, with rows 2 / distinct 1 printed: the "exactly one distinct" rule.
- **(N-absent):** a tmp ledger row naming a log that is absent from the tmp dir raises, naming it.
- **(N-digit):** a tmp copy of entry 62's log with the last digit of one lattice `sar10g` token changed gives exactly one (K) mismatch, naming that label. String equality resolves any digit, so this control has no ceiling issue.
- **(N-grid2):** a tmp copy with `:12407` duplicated raises.
- **(N-nogrid)** (entry 66's possible shape): a tmp copy with `:12407` deleted gives 51 fields, no `grid.*` names, and "grid: not run" printed.

**Printed.**
- (K)'s field count and its first and last names with their lines;
- the real tree's entry-66 status line — "(Pend)" today;
- the ledger's `(rows, distinct)` count.

**Windows** (bare and foreground, Bash `timeout: 660000`; `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 first; `FEM_EM_ANS_SKIP_PRIVATE=1` where the command carries it):
1. Code (and the guide sentence, if any). **Rule (i):** no page, test or guide edit after its last window.
2. **W1** `ANS-12-C58-W1`: C37's W1 command verbatim (`20261114T181121Z_ANS-11-C37-W1.log:12`: `-n 2`, `timeout -k 30 180`, the same five paths) with label `ANS-12-C58-W1` and raw log `ans12-c58-w1-raw.log`. Expect **84 passed per rank** (76 + 8 new), 0 failed / skipped, `[capture] rc=0`; count by token.
3. **Page** (the closing window): `scripts/testing/run_and_log.sh ANS-12-C58-page "./scripts/run_examples.sh -e ans:12 -n 1 -t 180"`. Expect 65 `[PASS]`, 0 `[FAIL]`, (R-cmp) 5051 B == 5051 B, and the (Pend) line. Then `git diff --stat` on ANS-12's `COMPARISON.md` and `ans:11`'s `metrics.json` / `COMPARISON.md` must be empty.
4. **Docrefs** `ANS-12-C58-docrefs`: `scripts/testing/run_and_log.sh ANS-12-C58-docrefs "timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py"` (item 983's form, `20261114T181259Z_ANS-11-C37-docrefs.log:12`). Expect `dead=0 guide=0`; exit 2 on the stale report is the precedent.

**Price / Tier / Class / Cost.**
- Standard, `-n 1` (the page; W1 at `-n 2`), complex build (the C6 import chain), no solve.
- **Class:** chain step (no solve). **Cost:** ≈ 10.8 measured (class median) / 40 predicted min.
- The windows sum to ≈ 35 s: C37's W1 24 s plus ≈ 3 s for the new nodes' tmp copies of a 1.6 MB log, the page 3 s, docrefs 2 s.

**Traps already paid for.**
- Items 871 / 930 / 931 / 983's lists: `re.search`, never `startswith`; either log suffix; derived values at ≤ 7 significant digits; open no `aed_results/`, `COMPARISON_private.md` or `docs/private/` file.
- **Do not edit the `ans:11` page or `test_ans11_record_page.py`.** Both are held by `xxl` entry 69 (premise 1). Do not edit `parse_c6_log`, `parse_c19_grid`, `_GRID_RE`, `custody_c19` or `CASES` either: call or read them only.
- **No top-level import added** to the `ans:12` page. The `ans:11` page's objects come through `_page()`.
- **The page import is collective** (`EX-193`): any test that loads it under `mpiexec` loads it on every rank.
- **The `GEO-38-C19` prefix:** match the chunk column exactly; never glob `*C19*`.
- **(Pend) is tree-independent:** never assert "absent" against the real ledger or log dir (items 905–907 → 918).
- **A pending window's reader parses its tracked log, never the host JSON** (`/workspace/logs/geo38/` is gitignored). (K) reads `metrics.json`, a tracked file.
- **A generator constant a landed test pins is never extended** (03:00 11-14): `D2_FIELDS_N` is new. Grep `D2_FIELDS_N` / `ENTRY66` across `tests/` first.
- Never write under `/workspace/logs/geo38` or `/workspace/logs/geo36`. Count pytest verdicts per rank, by token. Never delete a log, red or not.
- The leak hook: never `--no-verify`. On a refusal, move the log unchanged to `logs/leak-held/` and name the hook's line by number only.

**Scope.** The `ans:12` page's degree-2 reader and its "(Pend)" state, with (K) on entry 62. **Not claimed:**
- any 3 T degree-2 value (entry 66 has not run);
- ANS-12's `COMPARISON.md` degree-2 cells (they need the held `ans:11` page — after 11-28);
- the node (SPEC v2);
- the order rule on the 3 T column (entry 66's reading);
- any adjudication or AED number;
- the rod rows (entry 66's call 2; the C56 family).

**Negative result.**
- **(K) mismatch or (K-n) ≠ 53:** the `ans:12` route and C37's committed parse disagree on entry 62. That is a custody finding: stop, BLOCKED, and leave the known-issues entry to the review.
- **`parse_c6_log` raising on entry 62:** this contradicts C37's landed facts. Name the line, stop.
- **(R-cmp) red, or any landed node red:** something shared was edited. Stop, BLOCKED.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 C58 → DONE. `ANS-12` stays 🟡.
- Entry 66's reading becomes a slot step on its mechanical green, as the weekly wrote: the page parses the log the moment the ledger row names it.
- **Flag for the review:** filling ANS-12's `COMPARISON.md` degree-2 cells needs the `ans:11` page, which is held until `xxl` entry 69 runs. If entry 66 runs on 11-21 and entry 69 on 11-28, that half of entry 66's reading waits a week, or the weekly rules otherwise.

**Depends on:** nothing open. It is independent of item 988 (C56). It is never in the same slot as another item that edits the `ans:12` page or its test (none is queued).

**Outcome commit.** The logs, the two (or three) files, this file's Outcome paragraph, and the §9 line marked DONE with the local (CST) slot, in one commit through the hook.

**Outcome (12:00 slot, 2026-11-15, first try; GREEN).** The page gained `XXL_LEDGER`, `TRACKED_LOG_DIR`, `ENTRY66_CHUNK`, `ENTRY62_LOG_BASENAME`, `D2_FIELDS_N = 53`, `_ROW66_RE`, `entry66_scan` / `entry66_log`, `degree2_fields`, `fields_from_metrics`, `check_k_d2` and `render_d2`; no top-level import added, the `ans:11` page's parsers called unedited through `_page()`. The test gained 8 nodes (`test_c58_k`, `_pend`, `_stub66`, `_n_twin`, `_n_absent`, `_n_digit`, `_n_grid2`, `_nogrid`) with an own copy of C45's `_on_root`. Guide not edited.
- **W1** `20261115T181205Z_ANS-12-C58-W1.log` (`-n 2`, 27 s): **93 passed per rank**, 0 failed / skipped, `[capture] rc=0`. The brief's expected 84 counted C37's W1 (76); the other 9 are `test_ans11_record_page.py` nodes (34 → 43 per rank by node-token count, C37 W1 vs this log) that landed after C37's W1 (`79a8e180`, item 984, the only later commit to that file) — this slot did not edit that file. The `ans12` file has 17 nodes (9 + 8).
- **Page** `20261115T181338Z_ANS-12-C58-page.log` (`-n 1`, 2 s): **65 `[PASS]`, 0 `[FAIL]`**; (K-n) 53 vs 53, (K) 53 fields 0 mismatches, 0 unmatched metrics names, (K-src) `source_log` `20261114T080009Z_GEO-38-C19.log`, node 185.7629 (`:12160`); (R-cmp) live, 5051 B == 5051 B; (K) first field `s2.s11` (`:12379`), last `grid.sar10g` (`:12407`); the real-tree status line "(Pend) entry 66 pending — … (rows 0, distinct 0)". `git diff --stat` on ANS-12's `COMPARISON.md` and `ans:11`'s `metrics.json` / `COMPARISON.md` empty.
- **Docrefs** `20261115T181347Z_ANS-12-C58-docrefs.log` (2 s): `dead=0 guide=0`, exit 2 on the stale report (precedent).
- Controls, all green in W1: (N-twin) raises naming 2, the repeated row reads rows 2 / distinct 1; (N-absent) raises naming the stub base name; (N-digit) exactly one mismatch, `lattice.x-105mm dz-250mm.sar10g`; (N-grid2) raises on 2 grid lines; (N-nogrid) 51 fields, no `grid.*`, "grid: not run" rendered.
- Next: entry 66's reading becomes a slot step when the ledger names its row; ANS-12's `COMPARISON.md` degree-2 cells still wait on the `ans:11` page's hold (`xxl` entry 69).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 989 — `ANS-12` C58 — per `docs/planning/chunks/ANS-12.md` "§9 item 989". Edit only `examples/ansys_benchmarks/ans12_f_clinical_3t_box_loaded_128MHz/12_f_clinical_3t_box_loaded_128MHz.py`, `tests/unit/test_ans12_record_page.py` and, optionally, one digit-free sentence in the page's guide `.md`. Never edit the `ans:11` page or `test_ans11_record_page.py` (held by `xxl` entry 69), nor `parse_c6_log`, `parse_c19_grid`, `_GRID_RE`, `custody_c19` or `CASES`. ANS-12's `COMPARISON.md` and `ans:11`'s `metrics.json` / `COMPARISON.md` stay byte-identical. (K) compares string-exact against `ans:11`'s tracked `metrics.json`, never the host JSON; (Pend) runs on empty tmp roots. No solve, no mesh. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; container `timeout -k 30`; pytest `-s`. `pgrep -c python3` reads 0 first. Never edit a pattern, literal or band to rescue a red. On a red, park on `attempt/ANS-12-C58-*` and mark §9 item 989 BLOCKED in the same commit. Never delete a log. Never `--no-verify`. Do not spawn agents."
