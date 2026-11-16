# ANS-13 — the ASTM-box-loaded F-clinical birdcages' gel heating at 64 and 128 MHz against HFSS + IcepakFEA, one-way, per deposited watt

*History file for the §7 `ANS-13` row (`OPS-46` layout): the row holds the
state line and the Done-when; everything else lives here.*

## Commissioned — 2026-11-11 weekly planning review

Posed by `xl` entries 59 and 60's pre-registered decision rules ("every
assertion green ⇒ P33 lands as a record … and the 1.5 T `ANS` thermal case
on F-clinical (HFSS + Icepak, `ANS-9`'s pattern at clinical scale) is posed
by that weekly"; entry 60: "the 3 T thermal `ANS` case is posed beside the
3 T EM one"). Both windows ran with every executed assertion green —
entry 59 `20261110T080008Z_GEO-38-P33.log` (2026-11-10 02:00, 852 s,
41.70 GiB) and entry 60 `20261111T080008Z_GEO-38-P33-3T.log` (2026-11-11
02:00, 732 s, 32.75 GiB) — and are the record columns. The adjudication
columns are the degree-2 twins (`xl-pending.md` entries 69 / 70, PENDING on
§10 P87 and on the degree-2 roots), by the 2026-10-12 epitaph.

The SPEC is `examples/ansys_benchmarks/ans13_f_clinical_box_heating_64_128MHz/SPEC.md`
(v1). The EM half is `ANS-11` / `ANS-12` unchanged; the thermal half is
`MAT-15`'s route on the box (k 0.6, ρc 4.18e6, h 10 on six faces to 20 °C,
t\* 2821.5 s, CN N = 40), per watt deposited in the box; seven probes (c\*,
P0 … P5); readout bands `ANS-9`'s (5 % integral, 10 % pointwise with the
wall-probe ordering preserved).

**Steps (§10):** C53 — the `ans:13` record page (imports P61's / C31's
modules and C46's reader, parses entries 59 / 60's tracked logs; (X) / (D) /
(R) anchors; closing window `./run_examples.sh`); C54 — `COMPARISON.md`
through C36's script (the `ANS-11` / `ANS-12` rows byte-identical); P87 —
the degree-2 knob on P61's and C31's modules with the volume-mean θ(t\*)
print (entries 69 / 70's bodies; after entry 63 runs, the module hold);
entries 69 / 70 — the degree-2 columns.

**Not claimed at any outcome:** compliance, tissue, homogeneity, a dose per
watt (the per-watt `cem43` is the 20 °C baseline), a container wall, an
implant.

### §9 item 958 — 2026-11-12 (§10 C53, chain C — rule (4), the 11-11 weekly's step; lane 1 → chain step (no solve))

*Written by the 2026-11-12 03:00 daily review from the 11-11 weekly's C53 (`PROJECT_PLAN.md` §10 chain C), under §9 rule (4). It does not count against rule (6). Family `ANS-13` page, attempt 1 of 4. **Executor: implementer, in-slot.** This follows the `ans:12` precedent: item 930 (`ANS-12.md` § "§9 item 930", outcome DONE 2026-11-10 04:30) and item 833 before it were both run by the implementer. `example-runner` is for `EX-*` chunks. The 2026-11-05 03:00 ruling in `ANS-11.md` binds by analogy: a record reader, no solve, no mesh. The weekly's target binds: "entries 59 / 60's log tokens, C39's log strings and the SPEC v1 table".*

**PREDICTED RED (a SPEC premise, not a code defect). Read this before taking the item.** `SPEC.md:145`, 3 T column, gives the cw factor as `1.735230713e-01` with the cites (`:12294`, `:12304`). Entry 60's log prints the cw factor as `1.785856659e-01` at `20261111T080008Z_GEO-38-P33-3T.log:12293`. `:12304` is the "sources:" summary line, which carries no factor. C46 parsed and checked the log's value (`20261111T103024Z_GEO-38-C39-W1.log:356`). The SPEC token occurs in no tracked log. So (X) on that token and (R) on row `:145` are red by construction.
- **Step 0 of the slot:** if `SPEC.md:145` still reads `1.735230713e-01`, **do not take the item.** Journal "waits on the SPEC `:145` correction" and move on.
- The correction is the review's or the weekly's (the 11-11 weekly's precedent of correcting `ANS-12` SPEC v1's off-by-ones). The executor never edits the SPEC.
- Once the SPEC reads the log's token and cite, no `PREDICTED RED` remains.
- **Review note (2026-11-12 03:00 daily): corrected.** The review verified the drafter's reading (the token is in no log; `:12293` prints `1.785856659e-01`; C46 read it green at `…C39-W1.log:300`) and landed SPEC v1.0.1 errata in its own commit — `:145` now reads `1.785856659e-01 (`:12294`, `:12293`)`, and the two off-by-one citations (`:78` → `:12315`; `:99` → `_p37_fixture` `:155–157`) are fixed. The weekly ratifies or reverts. **Step 0 is therefore met; the item has no `PREDICTED RED`.** (R) regenerates the corrected table.

**Two readings of the step text against the files (execution detail; the target is unchanged).**
1. **"Found once" means matched once, to its row's parsed field** (the C48 reading, ratified 11-11). Raw multiplicity is printed. For example, C39's table is printed twice, from `:371` and from `:407` of the C39 log.
2. **"(R) byte-regeneration of the SPEC's record table" means `SPEC.md:140–148`:** the header, the separator and seven rows, Unicode minus and `…` included. The line citations outside that table are printed beside the parsed line, never asserted (the 930 disposition). Two are off:
   - the thermal table's "entry 59 `:12317`" for V / Robin area; the line is `:12315`;
   - row 2's "`test_geo36_rod_f_clinical_heating.py:443–445`, P37's points"; P37's points are `_p37_fixture` at `:155–157`, and `:443–445` is `_thermal_failures`.

   Both are SPEC transcription questions for the weekly.

**Code: new files only.** No SPEC edit, no edit to P61's, C31's or C46's module (P61's and C46's are held until entry 63 runs, 11-13 02:00, and P61's by entry 62 as well; **importing is fine**), no edit to the `ans:11` / `ans:12` pages, and no `src/` edit.
- **`examples/ansys_benchmarks/ans13_f_clinical_box_heating_64_128MHz/13_f_clinical_box_heating_64_128MHz.py`.** Pure functions, then `main()`, wrapped so that a raise calls `MPI.COMM_WORLD.Abort(1)` (the `ans:12` page's `:501`).
  - **Imports by name, never re-implemented.**
    - P61 (`tests/validation/test_geo38_f_clinical_box_heating.py`): `_probes`, `SENSE_TUPLE`, `T_STAR_BOX`, `RUNNING_WINDOW_S`, `HEAT_ENV`, `HEAT_SENSE_ENV`.
    - C31 (`tests/validation/test_geo38_f_clinical_3t_box_heating.py`) as `c31`: assert `c31.p61 is p61` (the 3 T column is P61's body unchanged under C21's patches, C31 docstring `:5–9`), and use `c31.C31_ENV`.
    - C46's reader (`tests/validation/test_geo38_c39_thermal_reader.py`): `find_log`, `parse_heat_log`, `check_block`, `side`, `ratios`, `_half_ulp_rel`, `_read_log`, `ENTRY59_GLOB` / `ENTRY60_GLOB`, `ENTRY59_NODE_PF` / `ENTRY60_NODE_PF`, `SENSES`.
    - Constants, by module name: `RHO_C`, `A_BOX` / `B_BOX` / `C_BOX` (`test_mat21_box_closed_form.py:31`, `:34`); `K_GEL`, `H_GEL`, `T_INF`, `N_MAIN` (`test_mat15_coil_driven_heating_phantom.py:82–88`); `T_STAR` (`test_mat21_box_fem.py:72`).
    - Never import a `test_*` name. Never call `_c39_body` (it writes C39's JSON at `both`) or `print_table`.
  - **Reads.**
    - Entry 59 is `find_log(ENTRY59_GLOB)` and entry 60 is `find_log(ENTRY60_GLOB)`. Each must return `"tracked"`.
    - Each log is parsed once by `parse_heat_log`: entry 59 is the block `:139` → gates `:12356`, entry 60 the block `:141` → gates `:12307`. Entry 60's cw probes sit after the cut at `:12298–12303` (§9 digest, item 955's note), and C46's per-line `re.search` reads them.
    - **New local regexes, applied only inside each block, for the lines C46 does not parse:**
      - the I32 source line's `V_sub` / `Robin area` (`:12315` / `:12265`);
      - the I32 (B) line `int_sub q` (ccw `:12316` / `:12266`, cw `:12334` / `:12284`);
      - the CN line's (E′) `stored … + surface …` and `(S2)` tokens (ccw `:12324` / `:12274`, cw `:12342` / `:12292`);
      - C6a's `normalisers:` line for `P_box/P_acc` and `SAR10g(c*) … per P_acc` (`:12312` / `:12262`).
    - Outside the blocks: the footer `Elapsed (s):` (`:12457` / `:12375`), the `memory.peak` bytes line (`:12452` / `:12370`), and each log's command line (`:12`) for the knob tokens.
    - **C39's log** `20261111T103024Z_GEO-38-C39-W1.log` is read with a local parser for the table rows under `C39 table:` (`:371`, 22 rows to `:393`). The second copy from `:407` is asserted identical, and the multiplicity is printed.
    - **A local SPEC parser** `parse_spec13(path)` works by regex, never by line number. It reads the record table (`:140–148`), the per-deposited-watt sentence (`:150–152`: 3.050e-02, 3.338e-02), the ratio sentence (`:153–156`: 4.0425 / 4.5397 with cites `:373` / `:384`, 4.1014, 3.6935), row 2's coordinates (`:99`), row 3's 9.737e-06 (`:100`) and the thermal table (`:78–86`).
  - **Prints.**
    - Every `[PASS]` / `[FAIL]` line.
    - Each probe and c\*, both senses, both field strengths: θ per P_acc watt, the factor, θ per deposited watt (the quotient, at `.6e`), and the I32 gel-twin token beside it.
    - The SPEC-line vs parsed-line table.
    - C33's `.6e` lever strings (`20261109T103557Z_GEO-38-C33-W1-r2.log:197–198`), printed only, beside the page's own (D-lever).
  - **Writes nothing.** (R-cmp) only reads.
  - **Rank safety:** C46's `find_log` / `parse_heat_log` are pure; call them on every rank, or through C46's `_on_root` on every rank, never under `if rank == 0` alone. If (R-cmp) is live, `c33._load_page()` is collective (`EX-193`, deviation 2). The page runs at `-n 1`.
- **Same-stem guide `.md`**, with the `ans:12` guide's headings:
  - the run line `./scripts/run_examples.sh -e ans:13 -n 1 -t 180`;
  - what is parsed, from which lines;
  - "degree-1 records at the degree-1 roots 193.5996 / 44.0525 pF; the adjudication columns are entries 69 / 70 (degree 2, PENDING on P87)";
  - per accepted watt vs per deposited watt (SPEC item 4);
  - links to `SPEC.md`, and to `COMPARISON.md` **only if item 959 has landed** (otherwise name it as C54's, unlinked).
  - It quotes closing-log digits only. Write `≈`, never a tilde.
- **Setup figure:** a `cmp`-identical copy of `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/figures/ans11_f_clinical_1p5t_box_loaded_64MHz_setup.png` (57 177 B, the 1.5 T coil with the box: the thermal domain is the box), saved as `figures/ans13_f_clinical_box_heating_64_128MHz_setup.png`. Disclose in the guide that the 3 T coil is `ans:12`'s figure. **`git add` it before W1** (item 930's W1 r1 went red 3 / 33 in `test_doc_reference_exit_codes.py` on an unstaged PNG).
- **`tests/unit/test_ans13_record_page.py`** (the `ans:12` test's loader shape, under its own `sys.modules` name) carries the gates and controls below.

**Gates (asserted in the script's `[PASS]` lines and in the unit test).**
- **(X) Extraction.**
  - **C46's `check_block(b, ENTRY5x_NODE_PF, list(p61._probes()), label)` on both blocks, every check green.** That covers (Node), (D-f), (P_gel), (D-scale), (Cem), (RA), (T\*) and (Probes), "the six probe names equal `p61._probes()` in order". C46 read every one of these green on the same two logs (`20261111T103024Z_GEO-38-C39-W1.log:336`, `:340–370`), so they are *asserted* (rule (e)).
  - SPEC row 2's seven coordinates equal `_probes()` values and the parsed hotspot point (`(-0.105, 0.0, -0.05)` / `(0.105, 0.0, 0.05)`) as floats. The SPEC prints `(−0.105, 0, −0.050)`, the 930 deviation.
  - **Custody:** entry 59's command line carries `p61.HEAT_ENV=1`, and entry 60's carries `c31.C31_ENV=1`. Both carry `p61.HEAT_SENSE_ENV=both` and their node.
  - **Every SPEC token of `:140–148` is matched to exactly one parsed field, by string, after normalisation** (strip a leading `+`, Unicode minus to `-`):
    - the factor, θ(c\*), P0 … P5, cw θ(c\*) and the cw factor;
    - the 1 g ball mean (the I32 ball token);
    - surface heat: the SPEC's 10-digit form equals `f"{float(tok):.9e}"` of the 13-digit (E′) token;
    - the 360 s running-mean maximum;
    - the header: 852 / 732 from the footers, 41.70 / 32.75 GiB as `f"{bytes / 2**30:.2f}"` (44777525248 / 35160018944), and the node strings;
    - every cell's line cite against the parsed line number.
  - **The thermal table (`:78–86`):**
    - k, ρc, h, T∞, t\* (`f"{T_STAR:.1f}"`), Δt (`T_STAR / N_MAIN`, and the log's `dt 70.537500`) and N equal the imported constants;
    - the box is 2 × (`B_BOX`, `C_BOX`, `A_BOX`) = 0.42 × 0.09 × 0.65;
    - V 0.02457 and Robin area 0.7386 are the `:12315` / `:12265` tokens. Their cite is printed (reading 2).
  - The one-match coverage line, "N of N", with none unmatched.
- **(D) Derived, at half a unit in the last printed digit.**
  - **(D-dep), 24 checks.** For each probe, both senses and both logs, θ per P_acc / factor equals the I32 gel-twin token. That token is θ at G·P = 1 W, i.e. per deposited watt (`:12318–12323`, `:12336–12341`; `:12268–12273`, `:12286–12291`). The bound is the Σ half-ulp of the three tokens. This is C46's (D-scale) in quotient form: equal to first order, and C46's tightest margin is 5.998e-10 ≤ 8.802e-10 (`…C39-W1.log:349`). *Asserted.*
  - **(D-c\*).** ccw θ(c\*) per P_acc / factor equals the SPEC's `3.050e-02` / `3.338e-02` (`:151`), at the Σ half-ulp, where the SPEC's 4-digit half-unit dominates (≈ 1.6e-4). Review arithmetic (unlogged): 3.04991e-02 and 3.33811e-02, rel 3.0e-5 and 3.3e-5.
    - **The cw c\* quotients have no printed comparator** in either log (I32 prints only the 1 g ball at c\*), so they are printed at `.6e`, never asserted. Review arithmetic: ≈ 3.0407e-02 and ≈ 3.2196e-02.
  - **(D-ratio).** C46's `ratios(side(b59), side(b60))` at `.6e` equals C39's 22 strings, by string. SPEC's `4.0425` / `4.5397` equal `f"{r:.4f}"` of the `:373` / `:384` ratios.
  - **(D-lever).** SPEC's `3.6935` equals `f"{fac60 / fac59:.4f}"` of the P61 ccw factors. SPEC's `4.1014` equals `f"{s60 / s59:.4f}"` of the normalisers' `SAR10g(c*) per P_acc` (2.274234e-03 / 9.327634e-03). Review arithmetic: 3.6934995 and 4.10144.
  - **(D-cal).** SPEC row 3's `9.737e-06` equals `f"{1 / (RHO_C * V):.3e}"`, with V the parsed `V_sub` token.
- **(R) Regeneration.** The page renders `SPEC.md:140–148` from the parsed tokens (Unicode minus, θ, `\*`, `…`, en-dash ranges), and the result equals those lines byte for byte.
- **(R-cmp), conditional.** If item 959 has landed (ANS-13's `COMPARISON.md` exists and the `ans:11` page has its `"ans13"` case), the page re-renders it through `c33._load_page()` and that case, and asserts byte identity. Otherwise it prints "absent — item 959 (C54)".
- `docrefs` `dead=0`.

**Controls (asserted).**
- **(N-digit-log):** in a `tmp_path` copy of entry 60's log, change the 4th significant digit of the ccw P3 per-P_acc token (`:12280`). (D-dep) goes red on that probe alone, and (X) / (R) go red. This is C46's (N-scale) shape (`test_c46_n_scale`).
- **(N-swap):** read entry 59's block as 3 T. (Node) is red, and (D-c\*) is red at ≈ 9.4 % against a ≈ 1.6e-4 bound.
- **(N-ratio):** a one-digit change to one ratio in a tmp copy of the C39 log makes (D-ratio) red.
- **(N-doc):** a one-digit change in a tmp SPEC copy's `:144` makes (R) red.

**Windows.** Bare, foreground, Bash `timeout: 660000`. `pgrep -c python3` in the container reads 0 first. `FEM_EM_ANS_SKIP_PRIVATE=1` on every pytest call.
0. **Pre-check:** `SPEC.md:145` carries the log's cw factor (see PREDICTED RED). Entry 63's queue file may still exist; that is fine (import only).
1. Census-pre: `scripts/testing/run_and_log.sh ANS-13-page-census-pre "./scripts/run_examples.sh --list"`. The expected count is 190 (`20261111T104016Z_EX-194-census-post.log`, 190 `-> examples/` lines), with ans 10. Read it; never assume it.
2. Docrefs-pre: `scripts/testing/run_and_log.sh ANS-13-page-docrefs-pre "timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py"`.
3. Code, guide and figure (staged).
4. **W1** `ANS-13-page-W1`: `scripts/testing/run_and_log.sh ANS-13-page-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans13_record_page.py tests/unit/test_doc_reference_exit_codes.py -v -s --tb=short > /workspace/logs/ans13-page-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans13-page-w1-raw.log; cat /workspace/logs/ans13-page-w1-raw.log; exit \$rc'"`. This is item 930's W1 form (`20261110T104239Z_ANS-12-page-W1-r2.log:12`; 33 passed, 23 s).
5. **Closing** `ANS-13-page`: `scripts/testing/run_and_log.sh ANS-13-page "./scripts/run_examples.sh -e ans:13 -n 1 -t 180"`. The `ans:12` closing took 2 s (`20261110T104306Z_ANS-12-page.log`). Nothing is edited after it, the guide included (rule (i)).
6. Census-post (191, ans 11), docrefs-post, and figures: `scripts/testing/run_and_log.sh ANS-13-page-figures "timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py"`.

**Tier.** Standard, `-n 1`, complex (the `ans:` group sources it; P61's import pulls in dolfinx), no solve, no mesh. **Class:** chain step (no solve), census +1. **Cost:** ≈ 5.1 measured (class median) / 45 predicted min. Comparable: item 930's windows ran 10:37:58 → 10:43:18 UTC (`…103758Z_ANS-12-page-census-pre.log` → `…104318Z_ANS-12-page-figures.log`), and this page parses two logs plus C39's. If W1 is not green by minute 45, park.

**Traps already paid for.**
- Pytest output is interleaved on log lines: entry 59 `:12343` and entry 60 `:12293` start with `PASSED…`, and entry 60's cw probes follow a test-id prefix at `:12298`. Use `re.search`, never `startswith`. Retention gzips logs: read either suffix through `_read_log`, and cite the base name.
- **Collective calls.** `_load_page()` on rank 0 alone hung `EX-193`'s first two W1s (`20261110T043240Z_EX-193-W1.log`).
- **The leak hook** (`MIN_SIGNIFICANT_DIGITS = 8`): derived prints use ≤ 7 significant digits, and raw log tokens are quoted as the `ans:11` page quotes them. On a refusal, move the log unchanged to `logs/leak-held/`, report it, and never use `--no-verify`.
- **Never call `_c39_body`.** Never write under `/workspace/logs/geo38` (C39's JSON). Tmp copies go to `tmp_path`.
- **The node is the log's token.** C31's docstring says 44.0615 (C16's stub node). Entry 60's node is 44.0525.
- **P61's "factor" is P_gel / P_acc.** With tag 403 counted as gel (`:12314`), that is the SPEC's P_box / P_acc. Print that identity; never re-derive the factor from the EM normaliser's 6-digit `P_box/P_acc`.
- Read no `aed_results/`, `COMPARISON_private.md` or `docs/private/`, and enter no AED value. The SPEC's status ticks are the weekly's.
- Never delete a log, a red one included.

**Scope.** The degree-1 record page on entries 59 / 60. **Not claimed:**
- the degree-2 columns (entries 69 / 70, P87);
- row 3's slope value: the logs print only its |relative miss| (S2);
- row 5 (P87 prints it);
- an adjudication or any AED number;
- a dose per watt, compliance, tissue, or a container wall.

**Negative result.**
- **(X) or (R) red on anything other than the pre-check's `:145`:** the parser or the log has moved, or a SPEC token is mistranscribed. Report the first mismatch, edit neither file, and mark BLOCKED. A SPEC mismatch goes to the weekly.
- **(D) red with (X) green:** a printed figure does not follow from its printed inputs. That is a finding on P61's / I32's writer. Park and mark BLOCKED; never retune.
- No other `PREDICTED RED:` line.

**Status it can move.**
- §10 C53 → DONE; census +1 (190 → 191, ans 10 → 11).
- **`ANS-13` Done-when (1)** (row `PROJECT_PLAN.md:3739`) needs this page **and** item 959's `COMPARISON.md` with the closing window. With 958 alone, (1) is "one `COMPARISON.md` (C54) away". With both, (1) is met.
- **The row's glyph is already 🟡** (set at SPEC v1, 11-11), so landing (1) moves the state line, not the glyph. ✅ is Done-when (2), the weekly's adjudication on entries 69 / 70. This brief reads the row's Done-when and matches it (rule (j)).

**Depends on:** the SPEC `:145` correction (pre-check). Nothing else: not item 959, and not entry 61 / 63 (it imports held modules, never edits them). If 959 has landed, (R-cmp) is live and the guide links `COMPARISON.md`.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 958 — `ANS-13` `ans:13` record page (C53) — per `docs/planning/chunks/ANS-13.md` "§9 item 958". First check `SPEC.md:145`; if it still reads 1.735230713e-01, stop and report. New files only (page, guide, unit test, figure copy); no SPEC, P61, C31, C46, `ans:11`, `ans:12` or `src/` edit. No solve, no mesh. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (2026-11-12, 04:30 slot; implementer): DONE, every gate green.** Step 0 met (SPEC v1.0.1 `:145` reads `1.785856659e-01`). Closing window `20261112T103734Z_ANS-13-page.log` (`-e ans:13 -n 1`, 3 s, exit 0): 165 `[PASS]`, 0 `[FAIL]`. (X): C46's `check_block` on both blocks (entry 59 `:139` → `:12356`, entry 60 `:141` → `:12307`), every check green; `c31.p61 is p61`; custody on both `:12` command lines; SPEC row 2's seven coordinates; the one-match coverage line "52 of 52, unmatched [], extra []" (tokens and cites); the thermal table's 11 checks (incl. the log's `dt 70.537500`, V / Robin area on all four `V_sub` lines). (D): (D-dep) 24 of 24; (D-c\*) 3.049908e-02 / 3.338108e-02 vs 3.050e-02 / 3.338e-02, rel 3.002e-05 / 3.248e-05 (≤ 1.639e-04 / 1.498e-04); cw c\* quotients 3.040714e-02 / 3.219548e-02 printed only; (D-ratio) 22 of 22 C39 strings, both copies (`:371`, `:407`) identical, 4.0425 / 4.5397 at the `:373` / `:384` cites; (D-lever) 3.6935 and 4.1014; (D-cal) 9.737e-06 on both `V_sub` tokens. (R) 1332 B byte-identical to `SPEC.md:140–148`. (R-cmp) printed "absent — item 959 (C54)". W1 `20261112T103702Z_ANS-13-page-W1.log`: 33 passed in 23.37 s (10 new nodes); controls red as predicted — (N-digit-log) (D-dep) red on 3 T ccw P3 alone (plus (D-ratio), which reads the same token), (X) `probes.p60[3]` and C46's (D-scale) on that probe, (R); (N-swap) both (Node) red, (D-c\*) rel 9.446e-02 / 8.631e-02; (N-ratio) (D-ratio) red; (N-doc) (R) red. Probe `20261112T103526Z_ANS-13-page-W0.log` (the script once, before the unit test; same 165 PASS). Census `…103044Z_ANS-13-page-census-pre.log` 190 (ans 10) → `…103738Z_ANS-13-page-census-post.log` 191 (ans 11); docrefs `dead=0` pre and post (exit 2, stale 122 unchanged, 190 → 191 guides); figures `…103742Z_ANS-13-page-figures.log` exit 0, `ans:13` ok. The SPEC v1.0.1 V / Robin cite `:12315` agrees with the parsed line (printed). **Deviations:** (R-cmp)'s live branch needs a hook in the `ans:11` page: once ANS-13's `COMPARISON.md` exists, this page loads that page through C33's `_load_page()` and requires a `render_ans13()` returning the file's text (no write). Item 959 should expose it, or (R-cmp) goes red. C33's lever strings are printed from their first copy (`:197–198`, of 2). §10 C53 is left to the weekly.

### §9 item 959 — 2026-11-12 (§10 C54, chain C — rule (4), the 11-11 weekly's step; lane 1 → chain step (no solve))

*Written by the 2026-11-12 03:00 daily review from the 11-11 weekly's C54 (`PROJECT_PLAN.md` §10 chain C), under §9 rule (4). It does not count against rule (6). Family C54, attempt 1 of 4. **Executor: implementer, in-slot** (the precedents are C36 / item 871 and C49 / item 931, `ANS-12.md` § "§9 item 931"). The weekly's target binds: "our degree-1 rows from entries 59 / 60 (rows 1, 2, 3, 4, 6 of the SPEC, both field strengths, cw at c\*), the degree-2 and AED columns blank, row 5 marked 'P87 prints it'; (K) the `ANS-11` / `ANS-12` rows regenerate byte-identical; `docrefs` 0 dead".*

**What "C36's comparison script" is, and its hold.** It is the `ans:11` page, `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/11_f_clinical_1p5t_box_loaded_64MHz.py`:
- `CASES` (`:653–697`), `render_comparison` (`:736–815`), `main_ans12` (`:818`) and `main(argv=None)` with `--case` (`:838–887`).
- C33's reader loads the page through `_load_page()`, and `mri:41` loads C33 (item 931's note).

**The page imports both of entry 62's modules at top level:** `test_geo38_f_clinical_box_loaded` (`:51`) and `test_geo36_rod_f_clinical_drive` (`:59`). Entry 62 is `xxl`, Saturday 11-14 02:00, `docs/testing/xxl-queue.d/20261114-GEO-38-C19.env`. By the import-closure hold (the 11-11 weekly's ratified ruling (d); §9 digest, "every module importing it … is held with it"), the page is held until entry 62 runs. The same holds for `tests/unit/test_ans11_record_page.py`, which loads it. **So this item is not takeable before entry 62's `chore(xxl)` commit is on `main` and its queue file is gone.** The step's "(now; serial on nothing)" predates that measurement. The hold is the binding ruling. Entry 62's window does not import the page, so the run itself is unaffected either way.

**Code: the `ans:11` page, `tests/unit/test_ans11_record_page.py`, and the new `examples/ansys_benchmarks/ans13_f_clinical_box_heating_64_128MHz/COMPARISON.md`. Nothing else.**
1. **A thermal case beside `CASES`, not inside it.** `render_comparison` / `build_metrics` read EM fields (`s2`, `y_in_ccw`, the lattice), so `"ans13"` goes into a separate `THERMAL_CASES` table, and `--case` choices become `sorted(CASES) + ["ans13"]`. The `"ans11"` / `"ans12"` paths stay byte-identical.
2. **`main_ans13(rank)`, with lazy imports.** The function imports C46's reader (`find_log`, `parse_heat_log`, `side`, `ENTRY59_GLOB` / `ENTRY60_GLOB`, `SENSES`) and P61's `_probes` **inside its body**, so the default import closure that C33 and `mri:41` load is unchanged.
   - It parses both entries, renders, writes ANS-13's `COMPARISON.md` only, and prints. No ANS-13 `metrics.json`, so `LISTED_EXAMPLE_ARTIFACTS` is unchanged (`.md` is not in `ARTIFACT_SUFFIXES`, item 931's outcome).
   - Local regexes, the same set as item 958's, for the lines C46 does not parse: (B), (E′) / (S2), `V_sub`.
   - **Never import the `ans:13` page** (it loads this page for (R-cmp): a cycle) **or C33's module** (the existing cycle).
3. **`render_comparison_thermal(blocks, extra)`.** It writes the title, the preamble (per deposited watt, SPEC item 4; degree 1 a record; the adjudication columns entries 69 / 70, PENDING on P87; AED never written here) and the table `| Row | Quantity | ours, degree 1, 1.5 T (entry 59, record) | ours, degree 1, 3 T (entry 60, record) | ours, degree 2 (entries 69 / 70, pending) | AED First Order + IcepakFEA |`. Each cell holds the raw token and its `:line`:
   - **Row 1:** the ccw (B) `int_sub q` token (`:12316` / `:12266`).
   - **Row 2:** at P0 … P5, the ccw I32 gel-twin θ(t\*) tokens (per deposited watt, `:12318–12323` / `:12268–12273`). At c\*, θ per P_acc / factor at `.6e`, with the 1 g ball token beside it. A **cw at c\*** row gives the cw quotient at `.6e`.
   - **Row 3:** the closed form `1 / (RHO_C · V_sub)` at `.4e`, with "(S2) |slope / closed form − 1|" and its token (`:12324` / `:12274`). The logs do not print the slope itself (`err_S2`, `test_geo36_rod_f_clinical_heating.py:350`), and the row says so.
   - **Row 4:** the ccw (E′) surface token (W per deposited watt).
   - **Row 5:** "P87 prints it".
   - **Row 6:** the ccw factor token (`:12325` / `:12275`), with the cw factor on the cw row.

   A provenance section names both logs (node, `-n 8`, 852 / 732 s, CN N = 40) and states what is not claimed.
4. **No degree-2 reader.** Entries 69 / 70 are PENDING, so those cells are blank.

**Gates (asserted in the unit test).**
- **(K) The `ANS-11` / `ANS-12` rows regenerate byte-identical**, against the files **as committed at the slot's start**: today `ans:11` `COMPARISON.md` is 4 566 B and `metrics.json` 6 025 B, and ANS-12's `COMPARISON.md` is 5 051 B. If C37 has landed first, re-read the sizes.
  - The default render of entry 52 and `build_metrics` equal the `ans:11` files.
  - `--case ans12` equals ANS-12's file.
  - After the closing windows, `git diff --stat` on all three is empty.
  - Every existing `test_ans11_record_page.py` node is green: C49's W1 read 67 passed, 1 skipped per rank (`20261110T104946Z_ANS-12-C49-W1.log`).
- **(X13).** Every degree-1 cell contains, for its row, the token string that entry 59's or 60's log prints, matched by `comparison_cells` (`:444`) or a thermal twin of it. Derived cells (the c\* quotients, the row-3 closed form) are recomputed from the parsed tokens and equal the cell's string.
- **(Probes).** Row 2's probe labels equal `list(p61._probes())` in order.
- **(Blank).** The degree-2 and AED cells are `""` on every row, and row 5's degree-1 cells read "P87 prints it".
- **(R13).** Re-rendering from the logs equals the committed ANS-13 `COMPARISON.md` byte for byte.
- `docrefs` `dead=0` (`COMPARISON.md` links `SPEC.md`).

**Controls (asserted).**
- **(N-case):** `--case ans13` rendered with entries 59 / 60 swapped fails (X13).
- **(N-digit):** a tmp log copy with the 4th digit of ccw P3's I32 token changed (`:12271`) fails (X13) on that cell alone.
- **(N-doc):** a one-digit change to a tmp copy of ANS-13's `COMPARISON.md` differs from the render.
- C36's and C49's existing controls stay as they are.

**Windows.** Bare, foreground, Bash `timeout: 660000`. `pgrep -c python3` reads 0 first.
0. **Pre-check:** entry 62's `chore(xxl)` commit is on `main`, and `docs/testing/xxl-queue.d/20261114-GEO-38-C19.env` is gone. Otherwise skip the item. No other slot is editing the `ans:11` page (C37).
1. Code.
2. **Wgen** `ANS-13-C54-gen`: `scripts/testing/run_and_log.sh ANS-13-C54-gen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/11_f_clinical_1p5t_box_loaded_64MHz.py --case ans13'"`. This is C49's Wgen form (`20261110T104926Z_ANS-12-C49-gen.log`, 2 s).
3. **W1** `ANS-13-C54-W1`: C49's W1 command (`20261110T104946Z_ANS-12-C49-W1.log:12`, `mpiexec -n 2`, `timeout -k 30 180`), labelled, with raw log `ans13-c54-w1-raw.log`. Add `tests/unit/test_ans13_record_page.py` to the paths if item 958 has landed. Leave C33's knob unset.
4. **Closing (K)** `ANS-13-C54-ans11-page`: `scripts/testing/run_and_log.sh ANS-13-C54-ans11-page "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"`. Then `git diff --stat` on `ans:11`'s `metrics.json` / `COMPARISON.md` must be empty.
5. **(K) 3 T** `ANS-13-C54-gen12`: Wgen's command with `--case ans12`. Then `git diff --stat` on ANS-12's `COMPARISON.md` must be empty. Then `scripts/testing/run_and_log.sh ANS-13-C54-ans12-page "./scripts/run_examples.sh -e ans:12 -n 1 -t 180"`, where (R-cmp) is live.
6. **Wreg** `ANS-13-C54-mri41`: `scripts/testing/run_and_log.sh ANS-13-C54-mri41 "./scripts/run_examples.sh -e mri:41 -n 2 -t 120"`. This is a regression only (C49's: 35 of 35, 2 s).
7. **If item 958 has landed:** `scripts/testing/run_and_log.sh ANS-13-C54-ans13-page "./scripts/run_examples.sh -e ans:13 -n 1 -t 180"`. (R-cmp) is then live and green. **This is the closing window of Done-when (1).**
8. **Docrefs** `ANS-13-C54-docrefs`: item 871's form, `dead=0`.

**Tier.** Standard, `-n 1` (`mri:41` at its own `-n 2`, W1 at `-n 2`), complex, no solve. **Class:** chain step (no solve). **Cost:** ≈ 5.1 measured (class median) / 35 predicted min. Comparable: item 931 ran Wgen → docrefs in 10:49:26 → 10:50:43 UTC; this item's windows sum to ≈ 2 min.

**Traps already paid for.**
- Item 871's / 931's list: `re.search`, never `startswith`; either log suffix; ≤ 7 significant digits on derived values; `FEM_EM_ANS_SKIP_PRIVATE=1`; a regenerated file is restored unless byte-identical.
- **No edit to `parse_c6_log`, `derived`, `parse_spec`, `CASES["ans11"]` or `CASES["ans12"]`:** C33 and `mri:41` read them.
- **No top-level import of C46 or P61 in this page.** It would widen the import closure that C33 / `mri:41` load. Import lazily inside `main_ans13`.
- **No import cycle:** this page imports neither C33 nor the `ans:13` page.
- C46's helpers are pure; any `_on_root` call runs on every rank. The page import is collective under `mpiexec` (`EX-193`).
- **C37 edits the same page** (entry 62's degree-2 rows). Never run both in one slot. (K) is against whatever is committed when the slot starts.
- The leak hook: never `--no-verify`. On a refusal, move the log unchanged to `logs/leak-held/`. No AED file is opened. Never delete a log.
- **No module or page is edited after its last window** (rule (i)). If the page changes after the `-e ans:11` window, re-run windows 4–7.

**Scope.** The comparison renderer's thermal case and ANS-13's `COMPARISON.md` with the two degree-1 columns. **Not claimed:**
- any degree-2 value or reader (entries 69 / 70, P87);
- row 3's slope value (only (S2) is printed);
- row 5 or row 7 as rows;
- an adjudication or any AED number;
- a SPEC edit.

**Negative result.**
- **(K) red:** the thermal case leaked into an EM render. Fix the default path, keep the first log, and re-run once (rule (i)).
- **(X13) red:** a parse or label error. Report it and stop, BLOCKED.
- **A C33 or `mri:41` regression red:** stop, BLOCKED, and leave the known-issues decision to the review.
- No `PREDICTED RED:` line. The SPEC `:145` error does not reach this item, because every cell comes from the logs.

**Status it can move.**
- §10 C54 → DONE.
- With item 958 landed, `ANS-13` Done-when (1) is met (`PROJECT_PLAN.md:3739`). The glyph stays 🟡, already set at SPEC v1, and the state line records (1). Otherwise `ANS-13` is "one page (C53) away", as the weekly wrote. This brief reads the row's Done-when and matches it (rule (j)).
- `ANS-11` / `ANS-12` stay 🟡.

**Depends on:** entry 62 having run (the hold). It does not depend on item 958: it reads entries 59 / 60 through C46, never through 958's page, and it may run before or after 958.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 959 — `ANS-13` C54, `COMPARISON.md` — per `docs/planning/chunks/ANS-13.md` "§9 item 959". First confirm entry 62's `chore(xxl)` commit is on `main` and its queue file is gone; otherwise stop. Edit only the `ans:11` page script and `tests/unit/test_ans11_record_page.py`, and write ANS-13's `COMPARISON.md` through the script's `--case ans13`; `ans:11`'s `metrics.json` / `COMPARISON.md` and ANS-12's `COMPARISON.md` must stay byte-identical. No solve. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (06:00 slot 2026-11-14): BLOCKED** (attempt 1 of 4, family C54; implementer, foreground). The pre-check was met: entry 62 is `3185acfd`, its queue file is gone and its Status is 0. The code was written: `THERMAL_CASES`, `render_ans13(log59=None, log60=None)` returning the text the `ans:13` page's `check_r_cmp` needs, `main_ans13` with lazy C46 / P61 / `RHO_C` imports, and ten `test_c54_*` nodes.
- **Wgen** `20261114T121455Z_ANS-13-C54-gen.log` (Status 0, 3 s): wrote ANS-13's `COMPARISON.md`, 3349 B.
- **W1** `20261114T121512Z_ANS-13-C54-W1.log` (`-n 2`, with `test_ans13_record_page.py`): 1 failed / 84 passed / 1 skipped per rank. All ten C54 nodes are green:
  - (X13) 13 rows, no misses (`:299`); (Probes) (`:301`); (Blank) (`:303`); (R13) 3349 B = 3349 B (`:309`);
  - controls (N-case) 24 misses (`:317`), (N-digit) on the 3 T P3 cell alone (`:321`), (N-doc) (`:325`);
  - C49's in-test (K) on 6025 B / 4566 B, and ANS-12 5051 B (`:271`, `:280`).
  
  **The red is the existing premise node `test_c19_absent`** (`:139`). Entry 62 wrote both degree-2 JSONs that `read_c19` looks for.
- **Window 4** `20261114T121606Z_ANS-13-C54-ans11-page.log`: every (K) `[PASS]` (`:70–87`). The page then rewrote `ans:11`'s `COMPARISON.md` with degree-2 cells, 16 lines changed (`:117`, `:120`); `metrics.json` is byte-identical. The file was restored by `git checkout`. Windows 5–8 did not run.

The red was not predicted, and it is not the brief's "(K) red: the thermal case leaked" case: it comes from the gitignored JSONs, which are C37's scope. So nothing was fixed. Filed in known-issues 2026-11-14. **Unblock:** C37 lands first, retiring or re-posing `test_c19_absent` and committing the degree-2 `ans:11` `COMPARISON.md`; 959 then re-poses with (K) against those files and takes the parked code. Deviations to carry into the re-pose:
- row 2's cw row carries the cw factor, and row 6 the ccw factor only;
- row 3's label is written "abs(…)" instead of using pipes;
- row 3 prints `9.7369e-06` at `.4e`, where the SPEC rounds to `9.737e-06`.

### Review ruling 2026-11-14 10:30 (item 959, C54 — BLOCKED on `test_c19_absent`, attempt 1 of 4) — re-posed as item 984 under rule (6)

**The red was the tree's, not the brief's or the slot's.** Entry 62 wrote the two host JSONs that C36's `read_c19()` looks for, and `main()` fed them into the default `ans:11` render. So `test_c19_absent` (a premise node) went red and `-e ans:11` rewrote `COMPARISON.md` from gitignored files. That is C37's scope. Item 983 (C37) moves the default render onto entry 62's **tracked** log, keeps the JSONs as custody only, and retires `test_c19_absent`. **Ruling:** C54 re-poses as item 984, **serial on 983 having landed in an earlier slot** (both edit the page; "never in the same slot"), with (K) against whatever `ans:11` / ANS-12 files are committed when 984's slot starts. Target, gates, controls and status-it-can-move are item 959's, unchanged. **The three disclosed deviations are accepted as consistent with the brief:** row 6 carries the ccw factor with the cw factor on the cw row (the brief's own wording); row 3's label says "abs(…)" because a `|` inside a Markdown table cell splits the cell; row 3 at `.4e` prints five significant digits (`9.7369e-06`), the brief's format. (X13) is against the log tokens, not the SPEC. The SPEC's 4-digit rounding is a SPEC-side print and needs no edit. Family C54: attempt 2 of 4. Provisional for the 11-14 weekly.

### §9 item 984 — 2026-11-14 (§10 C54, chain C — re-pose of item 959 under rule (6), attempt 2 of 4; chain step (no solve))

*Written by the 2026-11-14 10:30 daily review (re-pose of a blocked step, target kept). Provisional until the 11-14 weekly ratifies or reverts it. Family C54, attempt 2 of 4. **Executor: implementer, in-slot.** The ruling above and item 959's brief bind it.*

**Executes** §10 C54 exactly as item 959 states it (§ "§9 item 959" above: `THERMAL_CASES`, `render_comparison_thermal`, `main_ans13` with lazy imports, the six rows, the provenance section; gates (K), (X13), (Probes), (Blank), (R13), docrefs; controls (N-case), (N-digit), (N-doc); windows 0–8). **Only the start and (K)'s reference files change.**

**Pre-check (else skip at no cost):** item 983 (C37) is DONE on `main` in an **earlier** slot (its outcome commit and `ANS-11-C37-page-r` log present); `test_c37_custody` exists and `test_c19_absent` does not; known-issues 2026-11-14 is retired. Read and print the committed sizes of `ans:11`'s `metrics.json` / `COMPARISON.md` and ANS-12's `COMPARISON.md` at slot start — those, not 6 025 / 4 566 / 5 051 B, are (K)'s references.

**Start from the parked branch (`attempt/ANS-13-C54-20261114T1217Z`, `db34a506`, one commit on `04387052`).**
- **New files by path:** `git checkout attempt/ANS-13-C54-20261114T1217Z -- docs/testing/logs/20261114T121455Z_ANS-13-C54-gen.log docs/testing/logs/20261114T121512Z_ANS-13-C54-W1.log docs/testing/logs/20261114T121606Z_ANS-13-C54-ans11-page.log` (attempt 1's records; never delete a log). Append attempt 1's three `test-results.md` rows **by hand**, copied from the branch. Never path-checkout `test-results.md`.
- **The page script and `tests/unit/test_ans11_record_page.py`: never by path checkout.** C37 has edited both since `04387052`, so a path checkout would revert C37. Re-apply the parked hunks with the Edit tool from `git diff 04387052 db34a506 -- <file>` (read it with `git show`). They are additive: the thermal case, `render_ans13`, `main_ans13`, ten `test_c54_*` nodes, and one changed line, the `--case` choices. If a hunk's context no longer exists after C37, place it beside the same function and disclose the move.
- ANS-13's `COMPARISON.md` is regenerated by Wgen, never copied.

**Gates:** item 959's, with (K) against the slot-start files above. **Also green:** `test_c37_custody` and every C37 node (they are now landed gates in the same module).

**Windows:** item 959's 0–8, labels `ANS-13-C54b-*`. Window 4's `git diff --stat` on `ans:11`'s files must now be **empty**: C37 committed the tracked-log render, and that is the whole point of the re-pose. Window 7 (`-e ans:13`) is live, since item 958 landed (`f1ebddfc`), and it is the closing window of `ANS-13` Done-when (1).

**Tier.** Standard, `-n 1` (W1 and `mri:41` at `-n 2`), complex, no solve. **Class:** chain step (no solve). **Cost:** ≈ 10.8 measured / 30 predicted min (the code exists; attempt 1's windows summed ≈ 2 min).

**Traps already paid for.** Item 959's list verbatim, plus:
- **Path checkout of a file C37 edited reverts C37.** Take code by hunk, records by hand.
- A window that runs `-e ans:11` and leaves a `COMPARISON.md` diff after C37 means the default render still reads a host file. That is a finding: stop, BLOCKED. Never restore and continue.
- (K)'s byte sizes are re-read, never copied from item 959's text.

**Scope / negative result / status it can move:** item 959's, unchanged. §10 C54 → DONE; with item 958 landed, `ANS-13` Done-when (1) is met (glyph stays 🟡). The review deletes `attempt/ANS-13-C54-20261114T1217Z` once this lands. No `PREDICTED RED:` line.

**Depends on:** item 983 landed in an earlier slot. If it has not, skip at no cost.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 984 — `ANS-13` C54 re-posed — per `docs/planning/chunks/ANS-13.md` "§9 item 984", which imports item 959's brief. First confirm item 983 (C37) is DONE on `main` from an earlier slot and `test_c19_absent` is gone; otherwise stop. Take attempt 1's three logs by path from `attempt/ANS-13-C54-20261114T1217Z` and append their test-results rows by hand. Re-apply the parked page / test hunks with the Edit tool. Never path-checkout the page script, its test module or `test-results.md`. (K) is against the files committed at slot start. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. No solve. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (13:30 slot 2026-11-14): DONE** (attempt 2 of 4, family C54; implementer, in-slot, no executor spawned). Pre-check met: 983 is `6dcc2233` (12:00 slot), `test_c37_custody` present, `test_c19_absent` gone. (K)'s slot-start references: `ans:11` `COMPARISON.md` 9426 B, `metrics.json` 14192 B, ANS-12 `COMPARISON.md` 5051 B. Attempt 1's three logs taken by path from `attempt/ANS-13-C54-20261114T1217Z`, their three `test-results.md` rows appended by hand. The page and test hunks were re-applied by Edit. All context lines survived C37 unchanged, so nothing was moved, and the `-U0` added/removed lines of both files match `git diff 04387052 db34a506` byte for byte. C37 folded degree 2 into `build_metrics`, so `test_c54_k_cases_untouched`'s `render_comparison(m)` is still the default render.
- **Wgen** `20261114T193244Z_ANS-13-C54b-gen.log` (Status 0, 4 s): probes ok, wrote ANS-13 `COMPARISON.md` 3349 B, byte-identical to attempt 1's render.
- **W1** `20261114T193259Z_ANS-13-C54b-W1.log` (`-n 2`, 26 s, `[capture] rc=0`): 93 passed, 3 skipped per rank. (X13) 13 rows, misses `[]`; (Probes); (R13) 3349 B = 3349 B; (N-case) 24 misses; (N-digit) `:12271` misses the 3 T P3 cell alone. Every C37 node is green.
- **Window 4** `20261114T193334Z_ANS-13-C54b-ans11-page.log` (Status 0, 80 [PASS], no FAIL): `git diff --stat` on `ans:11`'s `metrics.json` / `COMPARISON.md` is **empty**.
- **Window 5** `20261114T193345Z_ANS-13-C54b-gen12.log` (Status 0, wrote 5051 B, `git diff --stat` empty) and `20261114T193355Z_ANS-13-C54b-ans12-page.log` (Status 0, 62 [PASS], (R-cmp) 5051 B = 5051 B).
- **Window 6** `20261114T193357Z_ANS-13-C54b-mri41.log`: 35 of 35 gates green.
- **Window 7 (closing, Done-when (1))** `20261114T193413Z_ANS-13-C54b-ans13-page.log`: Status 0, 166 [PASS], no FAIL, (R-cmp) live with 3349 B = 3349 B.
- **Docrefs** `20261114T193433Z_ANS-13-C54b-docrefs.log`: `RESULT: dead=0 guide=0 stale=121`, exit 2 (staleness only).

No file was edited after its windows. §10 C54 is DONE and `ANS-13` Done-when (1) is met (the glyph stays 🟡 on (2)). The review deletes `attempt/ANS-13-C54-20261114T1217Z`.
