# ANS-11 — the ASTM-box-loaded F-clinical-1.5T under the two-port drive at its Larmor root, 64 MHz

*History file for the §7 `ANS-11` row (`OPS-46`). Commissioned by the 2026-11-04 weekly (chain C's C7). The SPEC is `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/SPEC.md`. The adjudication column is the degree-2 one (C19, `xxl` entry 62); the degree-1 column (entry 52) is a record. AED values never enter this file.*

### Review ruling 2026-11-05 03:00 — what "the `ans:` page from entry 52's log, no run" means

The page is a **record reader**, on the `ans:4` item 417 precedent (`docs/planning/chunks/ANS-4.md:468–577`), which landed in one slot. It parses the tracked C6 log, and makes no solve and builds no mesh.
- **Why a reader:** C6 itself (868 s, 41 GiB, `-n 8`) cannot run in a slot. The gitignored `logs/geo38/c6-node-193.5996.json` is not reproducible from a fresh clone, so the tracked log is the source.
- **What it meets:** Done-when (1)'s degree-1 half only. The degree-2 rows wait on C19, so the glyph can move ⬜ → 🟡, not ✅.
- **The setup figure:** `write_setup_figure` needs a mesh (`src/fem_em_solver/post/setup_figure.py:180`), so the census carries this guide as `missing`. That is named for the figure lane.

### §9 item 833 — 2026-11-05 (lane 1, chain C step C7; `ANS-11` the `ans:11` record page, degree-1 column)

**Executes** the 11-04 weekly's lane-1 target "`ANS-11`'s runnable half" under the ruling above. Family `ANS-11` page, attempt 1 of 4.

**Code.**
- **New script** `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/11_f_clinical_1p5t_box_loaded_64MHz.py`, built from pure functions:
  - `_read_log`, copied from `ans:4` `:301–311`, with a `.log.gz` fallback.
  - `parse_c6_log(path)`. It `re.search`es `[GEO-38 C6a]` lines and keeps each token's **string** and 1-based line number:
    - `(Z) node` at `:12063` (V_s, Z_th);
    - `S2` at `:12281`;
    - `Y_in ccw` / `cw` at `:12285–12286`;
    - `mean |B1+|_ccw` at `:12288`;
    - the 18 `lattice` lines at `:12289–12306`;
    - `c*` at `:12307`;
    - `normalisers:` at `:12310`;
    - `wrote` at `:12311`;
    - `price:` at `:12312`.
  - It writes `metrics.json` (`json.dumps(indent=1, sort_keys=True)+"\n"`).
  - A `COMPARISON.md` writer for SPEC rows 1–6, with three columns:
    - ours, degree 1 (entry 52, a record, never adjudicated);
    - ours, degree 2 (C19, *pending*, blank);
    - AED First Order (blank).
- **Imports** from C6's module (`tests/validation/test_geo38_f_clinical_box_loaded.py`): `_lattice` (`:116`), `BRIEF_LATTICE_N` (`:51`), `C6A_LATTICE_DZ_M` (`:53`) and `_p_acc_spec` (`:376`). Never import a `test_*` name, and never re-implement a readout.
- **Private side:** no private writer and no `aed_results/` read. Honour `FEM_EM_ANS_SKIP_PRIVATE`.
- **Same-stem guide `.md`:**
  - the run line `./scripts/run_examples.sh -e ans:11 -n 1 -t 180`;
  - what is parsed, and from which lines;
  - "a degree-1 record; the adjudication column is C19";
  - links to SPEC.md and COMPARISON.md.
- **Other files:**
  - Declare `metrics.json` in `LISTED_EXAMPLE_ARTIFACTS` (`tests/unit/test_doc_reference_exit_codes.py:473`) in the same commit.
  - New `tests/unit/test_ans11_record_page.py`.
  - No SPEC edit; the SPEC's status ticks are the weekly's.

**Windows.** Bare, foreground, Bash timeout 660000.
1. `scripts/testing/run_and_log.sh ANS-11-page-list-pre "./scripts/run_examples.sh --list"` (ans 8).
2. Code.
3. W1: `scripts/testing/run_and_log.sh ANS-11-page-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans11_record_page.py tests/unit/test_doc_reference_exit_codes.py -v -s --tb=short'"`
4. Closing window: `scripts/testing/run_and_log.sh ANS-11-page "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"`
5. `ANS-11-page-list-post` (ans 9, census 175 → 176).
6. `ANS-11-page-docrefs` (`timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py`, dead = 0).
7. `ANS-11-page-figures` (`timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py`). It reads `missing` +1, exit 2; name it.

**Anchor (asserted in the script and in W1).**
- **(X) Extraction.**
  - The 18 lattice lines' labels equal `_lattice()`'s, in order.
  - Every SPEC-quoted value equals its parsed token: S, Y_in, the means, c\*, and P_box / P_src / P_acc (`SPEC.md:117–123`).
  - The 36 six-digit lattice entries (`:129–131`) equal `f"{float(tok):.6e}"`.
  - SPEC values are parsed from SPEC.md, never retyped.
- **(D) Derived quantities** from the parsed tokens, each at half a unit in its last printed digit:
  - `_p_acc_spec(S2, (1, −1j), V_s, 50.0)` against SPEC's P_acc (rel ≤ 1e-8);
  - SAR10g(c\*)/P_acc against 2.274234e-03 (≤ 5e-10);
  - P_box/P_acc against 5.285093e-02 (≤ 5e-9);
  - the 18-point lattice mean against the SPEC's 1.49792e-08 (≤ 5e-14);
  - the argmax of the 18 SAR values equals the c\* label, and the c\* SAR string equals that lattice line's string;
  - the S₁₂ string equals the S₂₁ string.
- **(R)** The writer, run on the committed `metrics.json`, reproduces the committed `COMPARISON.md` byte for byte.

**Control (asserted).**
- **(N-digit-log):** in a tmp copy of the log, change one digit of the c\* lattice line's SAR. (D)'s c\* string equality fails, and so does (X).
- **(N-drop):** delete one lattice line. The parser raises (17 ≠ `BRIEF_LATTICE_N`).
- **(N-doc):** a one-digit change to a tmp `COMPARISON.md` differs from the render.
- **(N-sense):** the cw weights (1, +1j) miss P_acc. The size is predicted from hand arithmetic, so it is printed, not asserted (rule (e)).

**Tier.** Smoke / standard, `-n 1`, complex build (the C6 import pulls in dolfinx). **Class:** example. **Cost:** ≈ 7.3 measured / 50 predicted min. This is item 417's footprint plus a guide. If W1 is not green by minute 45, park.

**Traps already paid for.**
- Pytest output is interleaved on log lines (e.g. `:12321 …SKIPPEDPASSED`). Use `re.search`, never `startswith`.
- Retention gzips logs: read either suffix, and cite the base name.
- The private-leak commit hook (`MIN_SIGNIFICANT_DIGITS = 8`) can flag a derived digit string, as it did item 417's class means. If it does, print fewer digits and move the first log, unchanged, to `logs/leak-held/`.
- Declare the artifact in `LISTED_EXAMPLE_ARTIFACTS`.
- Never run C6, never hand-edit `COMPARISON.md`, and enter no AED value.
- The 400-centroid mean is a second record beside the lattice mean, not the comparison row (the 11-04 weekly's correction (i)).

**Scope.** Degree-1 record rows only. **Not claimed:** the degree-2 rows, an adjudication or any AED number, the node (SPEC v2's), or a setup figure.

**Negative result.**
- (X) red: the parser is wrong, or the log moved. Report the first mismatch and stop.
- (D) red with (X) green: a printed figure does not follow from its printed inputs (P_acc's convention, say). That is a finding for the review: park, BLOCKED, never retune.
- A mismatch between the SPEC and the log is a SPEC transcription error for the weekly. Do not edit the SPEC.
- No `PREDICTED RED:` line.

**Status it can move.**
- `ANS-11` ⬜ → 🟡: Done-when (1) is met for the degree-1 column. The degree-2 rows land with entry 62.
- Census 175 → 176 (ans 8 → 9).

**Depends on:** nothing. It does not wait for entry 56; the node changes only SPEC v2.

**Outcome (2026-11-05, 16:30 slot, implementer): DONE.** The page parses `20261103T080009Z_GEO-38-C6.log` at the cited lines; (X) green (18 labels, 20 SPEC values, 36 / 36 six-digit entries), (D) green (P_acc rel 1.255e-10; SAR10g(c\*)/P_acc 3.2e-10, P_box/P_acc 3.7e-9, lattice mean 2.7e-14 inside half a printed unit; argmax = c\*), (R) 4 566 B identical; (N-digit-log), (N-drop), (N-doc) caught; (N-sense) cw/ccw P_acc ratio 1.0053, printed. Closing window `20261105T223815Z_ANS-11-page.log` (30 PASS, 3 s); W1 `…223748Z_ANS-11-page-W1d.log` (28 passed, 2 failed = the filed dead=2 SPEC entry, which turns two `test_doc_reference_exit_codes.py` nodes red — added to that entry). Census 176 → 177 (ans 8 → 9; the brief's 175 was one behind); figures `missing=1` is this page (no mesh, so no setup figure). `ANS-11` ⬜ → 🟡.

### §9 item 845 — 2026-11-05 (ops / docs; the "Resolves with" of known-issues 2026-11-05, the `ANS-11` SPEC's two dead references)

**Why first.** The SPEC's bare basenames `build.py` / `box.py` make the docrefs checker read `dead=2` (exit 1) and turn two `tests/unit/test_doc_reference_exit_codes.py` nodes red (`20261105T223748Z_ANS-11-page-W1d.log:919–921`), so `main` has two known red test nodes and every example closing since 11-05 has had to disclose `dead=2` (`EX-182`, item 833, batches 4–5). **Identified by this review** (read, not executed): the eleven citations `build.py:46` … `build.py:59` are `tests/validation/test_geo38_f_clinical_build.py` (line 59 `shield_radius=0.3725`, 58 `air_padding`, 55 `port_box_size`, 46–51 `ring_radius` … `ring_gap_length`) and the two `box.py:51` / `box.py:60–62` are `tests/validation/test_geo38_f_clinical_box.py` (line 51 `TABLE_OFFSET_M = 0.0`, 60–62 `phantom_box` `size_m`). The SPEC already cites `tests/validation/test_geo38_f_clinical_build.py:44–64` in full at `:44`.

**Edit.** In `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/SPEC.md`, on lines 51–60 only, replace each `` `build.py:N` `` with `` `test_geo38_f_clinical_build.py:N` `` and each `` `box.py:N` `` with `` `test_geo38_f_clinical_box.py:N` `` (the basename the checker resolves; same line numbers). **Before rewriting, read each cited line in its file and confirm it holds the value the SPEC sentence quotes**; if any does not, stop and journal (a wrong citation is the weekly's to correct, not a rename). No number, no SPEC version and no line count changes (the page's parser reads `SPEC.md:117–131`; assert `wc -l` is unchanged). No other file edited. The known-issues entry is removed in the same commit.

**Anchor.** The checker's own contract (`OPS-19`): `RESULT: dead=0` and exit ≠ 1; the two named nodes green; `test_ans11_record_page.py` still green with its (X) / (D) / (R) gates (the page reads values, not citations).

**Control.** The same checker on the pre-edit tree read `dead=2` naming exactly these two basenames (`20261105T193155Z_EX-179-b5-docrefs.log`, the `RESULT:` line `:156`); the post-edit log names neither (`grep -c 'build.py\|box.py: no such file'` 0). The stale count is information and is not compared.

**Windows.** Bare, foreground, Bash timeout 660000 ms.
- `scripts/testing/run_and_log.sh ANS-11-spec-docrefs "timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py"`
- `scripts/testing/run_and_log.sh ANS-11-spec-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans11_record_page.py tests/unit/test_doc_reference_exit_codes.py -v -s --tb=short'"` — item 833's W1d form (`…W1d.log:12`); expect 30 passed, 0 failed.

**Tier.** Smoke / standard, `-n 1`, no solve. Predicted ≈ 25 s of windows. **Class:** ops / docs. **Cost:** ≈ 3.4 measured / 10 predicted min.

**Traps already paid for.** A harness call is the whole Bash command; `FEM_EM_ANS_SKIP_PRIVATE=1` on any `ans:` run; no `docs/private/` or `aed_results/` read; quote no AED number; count verdicts per rank by token; a regenerated `metrics.json` / `COMPARISON.md` is restored unless byte-identical.

**Scope.** Two citations' spelling. **Not claimed:** any SPEC content, SPEC v2, the `ans:11` setup figure (item 855).

**Negative result.** `dead` still > 0: print the checker's FAIL lines, stop, BLOCKED (the checker's basename rule is then the finding). A cited line that does not hold its quoted value: journal it for the weekly, edit nothing. No `PREDICTED RED:` line.

**Status it can move.** Known-issues 2026-11-05 (the SPEC's dead references) retired; §9's red-node bullet back to "no known red test node".

**Depends on:** nothing.

**Outcome (2026-11-05, 19:30 slot, implementer): DONE.** Each cited line read first and found to hold its quoted value (build `:46`–`:51`, `:55`, `:58`, `:59`; box `:51`, `:60–62`); eleven `` `build.py:N` `` → `` `test_geo38_f_clinical_build.py:N` `` and two `` `box.py:N` `` → `` `test_geo38_f_clinical_box.py:N` ``, SPEC lines 51–60 only, `wc -l` 153 before and after. Docrefs `RESULT: dead=0 guide=0 stale=116 … exit=2` (`20261106T013105Z_ANS-11-spec-docrefs.log:156`; the two basenames 0 times). W1 `20261106T013115Z_ANS-11-spec-W1.log`: 30 passed, 0 failed, 24 s, both named nodes PASSED, `test_ans11_record_page.py` green. Known-issues 2026-11-05 entry removed.

### §9 item 871 — 2026-11-08 (§10 C36, chain C — ANS-11's comparison script; lane 1 → chain step (no solve), standard, -n 1)

*Written by the 2026-11-08 03:00 daily review from the 11-07 weekly's enumerated C36 (`PROJECT_PLAN.md` §10 chain C, "C36 (now; lane 1 → chain step, standard, `-n 1`, no solve; serial on nothing)"), under §9 rule (4). Family C36, attempt 1 of 4. The 2026-11-05 03:00 ruling above binds: the page is a record reader, makes no solve and builds no mesh.*

**The weekly's text (binds).** "`COMPARISON.md`'s rows generated from a reader over entry 52's tracked log (the degree-1 record column, as the `ans:11` page already parses it) and over the JSON names entry 62 will write (`h9p-node-185.7629-d2.json` and the grid / rod JSONs; blank until they exist), the AED columns blank; (K) the degree-1 rows equal the SPEC's table by string; (D) the SPEC v2 node 185.7629 pF and R_p ‖ C figures re-derived from the module constants to the printed digit; `docrefs` 0 dead. Landing ⇒ C37 is one slot after entry 62."

**The JSON names, read against the code (a correction of spelling, not of target).** No module writes `h9p-node-185.7629-d2.json` (that is `WF-10`'s F-human family, `test_wf10_f_human_tuned_b1_sar.py:165`). Entry 62's two calls (`docs/testing/xl-pending.md:2453`) write:
- call 1, C6's module at degree 2: `/workspace/logs/geo38/c6-node-185.7629-d2.json` (`test_geo38_f_clinical_box_loaded.py:739–741`, `C6_JSON_DIR` `:79`). The box-wide grid is **a `"grid"` key inside that file** (`:756–764`, present only with `FEM_EM_GEO38_C6_GRID=1`, which call 1 sets), not a separate JSON. C6's file has no `"degree"` key (`:740`, "no new key"); the `-d2` suffix is the only degree mark.
- call 2, I36b's module at degree 2: `/workspace/logs/geo36/i36-node-185.7629-d2.json` (`_i36_json_path`, `test_geo36_rod_f_clinical_drive.py:219–224`; `I36_JSON_DIR` `:66`), which carries `"degree": 2` (`:610`).

These are the "grid / rod JSONs". Both are gitignored, and neither exists before 2026-11-14.

**Code: `examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/11_f_clinical_1p5t_box_loaded_64MHz.py` and `tests/unit/test_ans11_record_page.py` only.** Pure functions, added beside the existing ones. No `src/` edit, no SPEC edit, no edit to C6's, I36b's, C18's or `WF-10`'s module.
1. **`read_c19(c_pf="185.7629", c6_dir=None, i36_dir=None)`.** The paths are `os.path.join(c6_dir or C6_JSON_DIR, f"c6-node-{c_pf}-d2.json")` and `i36._i36_json_path(float(c_pf), 2)`, the latter re-rooted on `i36_dir` when one is given. Import `C6_JSON_DIR` and `K8_JSON_KEYS` (`test_geo38_f_clinical_box_loaded.py:79`, `:1231`) and `_i36_json_path` / `I36_JSON_DIR`. It returns `{"c6": rec | None, "i36": rec | None, "paths": [...]}`. An absent file gives `None` and no error. A present file is validated:
   - C6: `set(rec) == K8_JSON_KEYS` or `K8_JSON_KEYS | {"grid"}`, and `f"{rec['c_pf']:.4f}" == c_pf`;
   - I36: `rec["degree"] == 2`, and `f"{rec['c_pf']:.4f}" == c_pf`.
   Anything else raises `ValueError` naming the file and the key.
2. **`degree2_cells(c19)`.** Returns `{row label: cell string}` for the degree-2 column, with every cell `""` when `c19["c6"]` is `None`. The mapping, from the C6 JSON keys `:742–755`:
   - S₁₁, S₁₂ = S₂₁, S₂₂ ← `s2`;
   - Y_in ccw P33 / P41 ← `y_in["ccw"]`;
   - the 18-point lattice mean and its value per √P_acc ← the mean of `lattice[*]["b_plus_ccw"]`, and `p_acc`;
   - c\* ← `star_label`; SAR10g(c\*) and its value per P_acc ← `sar_star`, `p_acc`; the cw ball ← `sar_star_cw`;
   - the box-wide 10 g maximum ← `grid["star_label"]` / `grid["sar_star"]` when `"grid"` is present;
   - the 400-centroid mean and its value per √P_acc ← `mean_b_plus["ccw"]`;
   - P_box / P_src / P_acc; gel fraction P_box / P_acc;
   - sense ← `mean_b_plus["ccw"] / mean_b_minus["ccw"]`.
   Raw values use the C6 module's own print formats (`.10e` for S as at `:12281`, `.9e` elsewhere). Derived values use the page's existing derived formats (`.6e` / `.5e` / `.4e`, ≤ 7 significant digits, the leak hook's floor), so that C37's degree-1 → degree-2 move is a like-for-like string per row. The copper / ESR fractions row stays `""` (C37's to decide). **The I36 record maps to no row:** the SPEC has no rod row ("Out of scope: no implant"). It is read and validated only.
3. **`render_comparison(m, c19=None)`** takes its degree-2 cells from `degree2_cells(c19)`. With `c19=None`, or with both files absent, the output is **byte-identical** to today's committed `COMPARISON.md` (4 566 B). `build_metrics` is unchanged when absent (`"degree2_c19": "pending (xxl entry 62)"`). `main()` calls `read_c19()` and prints one line naming both paths and `present` / `absent`. It writes nothing new.
4. **(K) `check_k(spec, comparison_text)`.** It parses the degree-1 cell (third column) of each COMPARISON row by its Quantity label. The normalisation N is: drop every `(…)` parenthetical, delete spaces, map `−` → `-`, strip one leading `+`. It asserts that N(SPEC token) occurs in N(cell) for every token `parse_spec` returns:
   - `s2` (s11, s12, s22) and `y_in_ccw` (p33, p41);
   - `means.b1p_ccw`, `means.sense` and `lattice_mean`;
   - `c_star` (sar10g, sar10g_cw);
   - `powers` (p_box, p_src, p_acc);
   - `per_p_acc` (all three).
   It also asserts that the SPEC's 36 six-digit lattice entries equal `f"{float(cell):.6e}"` of COMPARISON's lattice table, matched by label. SPEC tokens that COMPARISON does not carry (`means.b1m_ccw`, `means.cv`, and `c_star`'s x / z / r / az in the SPEC's form) are **printed as "not carried"**, and their count is asserted equal to that list of 6, so no omission can be silent. The SPEC table sits at `SPEC.md:120–126` and `:132–134` since v2. `parse_spec` reads by regex, not by line, but the script's docstrings still cite v1's `:117–123` / `:129–131`. Correct those two citations in the same edit, disclosed.
5. **(D) `node_figures(c_f, frequency_hz=FREQUENCY_HZ_C4, r_ohm=ESR_OHM, z0=50.0)`.** The constants are `FREQUENCY_HZ_C4` (`test_geo38_f_clinical_esr_root.py:67`) and `ESR_OHM` (`test_wf10_f_human_copper_esr_root.py:53`). Z_RC comes from that module's `_z_rc` (`:188–191`) and Z_th from its `_thevenin_feed_impedance` (`:144`); never re-implement either. With ω = 2π·f, it returns:
   - Q_C = 1/(ωCr);
   - R_p by the SPEC's stated high-Q form 1/(ω²C²r) (`SPEC.md:75`);
   - C_p = C/(1 + Q_C⁻²);
   - Z_RC and Z_th.
   `parse_spec_node(path)` parses, by regex and never by retyping:
   - the node string at `SPEC.md:14` and its repeat at `:75`;
   - `:75`'s v2 figures (R_p, the C beside it, Q_C, C_p, Z_RC);
   - `:75`'s v1 parenthetical (R_p, C, Q_C, C_p, Z_RC);
   - `:77`'s v2 Z_th and the v1 Z_th beside it.
   Compare each derived value to the SPEC's string at that string's own printed decimals.

**Gates (asserted, in the unit test and in the script's `[PASS]` lines).**
- **(R) unchanged.** The committed `metrics.json` (6 025 B) and `COMPARISON.md` (4 566 B) regenerate byte for byte. `git diff --stat` after the page window is empty on both files.
- **(C19-absent).** On this tree both entry-62 paths are absent. `read_c19()` returns `None` / `None`, and every degree-2 cell is `""`.
- **(C19-stub)**, on `tmp_path` with rank 0 writing and a bcast. A C6 stub has exactly `K8_JSON_KEYS ∪ {"grid"}` (values are arbitrary floats; the two-sense `y_in` / `mean_b_plus` shape follows `:745–748`). An I36 stub has `{"c_pf": 185.7629, "degree": 2}`. On these stubs:
  - each mapped degree-2 cell is non-empty;
  - every degree-1 cell is unchanged;
  - (K) is still green on the stub render, which shows that (K) reads only the degree-1 column.
- **(K)** green on the committed `COMPARISON.md`, with the six "not carried" tokens printed.
- **(D-node)** The SPEC's node string (`:14`, `:75`) equals `f"{ROOT_D2_PF:.4f}"` = `185.7629` (`test_geo38_c18_degree2_c_mode.py:100`, entry 56's crossing 185.762904 pF).
- **(D-v1), the anchor that ties the page's arithmetic to the module's own print.** At C = 193.5996 pF (the C6 log's `node` token, `20261103T080009Z_GEO-38-C6.log:12063`):
  - `f"{z_rc:.9e}"` and `f"{z_th:.9e}"` equal that line's `Z_RC 2.000000000e-01-1.284504713e+01j` and `Z_th 3.259458821e+00-1.195984969e+01j` strings;
  - R_p 825.0, Q_C 64.2, C_p 193.5527 and Z_RC −12.845047127j equal the SPEC's v1 parenthetical (`:75`).
- **(D-v2).** R_p, Q_C, C_p, Z_RC and Z_th equal the SPEC's v2 strings (`:75`, `:77`: 896.1, 66.9, 185.7215, −13.386935347j, 3.505604063 − 12.39875443j), derived **at C = `ROOT_D2_PF` = 185.762904 pF**.
- **Printed, not asserted (rule (e); review arithmetic, unlogged):** the same figures at the **node** C = 185.7629 pF, the value entered in HFSS and in entry 62's `FEM_EM_GEO36_S8_DRIVE_C_PF`. Predicted: R_p 896.1 and Q_C 66.9 unchanged, but C_p **185.7214**, Z_RC −13.386935636j and Z_th 3.505604196 − 12.39875466j, so they move in the last printed digit(s). Printed beside them: R_p by the binding ruling's exact form r(1 + Q²), predicted 896.3 at either C, where the SPEC prints the high-Q form.

**Controls (asserted).**
- **(N-K)** One digit of S₁₁ changed in a tmp copy of SPEC.md: (K) goes red naming S₁₁, and (X) goes red too.
- **(N-D)** The v2 figures derived at C = 193.5996 pF miss SPEC v2's R_p string (825.0 ≠ 896.1).
- **(N-stub-key)** A C6 stub with one `K8_JSON_KEYS` key removed makes `read_c19` raise `ValueError` naming that key.
- **(N-stub-node)** A C6 stub whose `c_pf` is 193.5996, written under the 185.7629 name, raises.
- **(N-stub-deg)** An I36 stub with `"degree": 1` raises.
- **(N-digit-log), (N-drop), (N-doc)** stay as they are.

**Windows.** Bare, foreground, Bash timeout 660000 ms. No compute runs beside them; `pgrep -c python3` in the container reads 0 first.
- **W1:** `scripts/testing/run_and_log.sh ANS-11-C36-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans11_record_page.py tests/unit/test_doc_reference_exit_codes.py -v -s --tb=short > /workspace/logs/ans11-c36-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans11-c36-w1-raw.log; cat /workspace/logs/ans11-c36-w1-raw.log; exit \$rc'"`. This is item 845's W1 (`20261106T013115Z_ANS-11-spec-W1.log:12`; 30 passed, 24 s, `:905`, `:909`) inside the durable capture of `20261015T171259Z_GEO-38-s1a-W2-r2-a2.log:12`, escaped as `GEO-38.md:233` spells it. Expect 30 passed plus the new nodes, 0 failed.
- **Page:** `scripts/testing/run_and_log.sh ANS-11-C36-page "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"` (3 s, `20261105T223815Z_ANS-11-page.log:78`). Every `[PASS]` line is followed by the new `read_c19` line (`absent`). Run `git diff --stat` on `metrics.json` / `COMPARISON.md` afterwards; it must be empty.
- **Docrefs:** `scripts/testing/run_and_log.sh ANS-11-C36-docrefs "timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py"`, which must read `dead=0` (item 845's `20261106T013105Z_ANS-11-spec-docrefs.log:156`).

**Tier.** Standard, `-n 1`, complex (the C6 import pulls in dolfinx), no solve. **Class:** chain step (no solve). **Cost:** ≈ 4.8 measured / 30 predicted min. The windows sum to ≈ 35 s; the hour goes on code.

**Traps already paid for.**
- Pytest output is interleaved: use `re.search`, never `startswith` (item 833).
- Retention gzips logs, so read either suffix and cite the base name.
- The private-leak hook (`MIN_SIGNIFICANT_DIGITS = 8`) can flag a derived digit string. Print derived values at ≤ 7 digits. On a hit, move the log unchanged to `logs/leak-held/`, report, and **never `--no-verify`**.
- `FEM_EM_ANS_SKIP_PRIVATE=1` on every `ans:` run. Read no `aed_results/`, `COMPARISON_private.md` or `docs/private/`, and enter no AED value.
- A regenerated `metrics.json` / `COMPARISON.md` is restored unless it is byte-identical.
- Stubs live on `tmp_path`: rank 0 writes, then a bcast. **Never write under `/workspace/logs/geo38` or `/workspace/logs/geo36`**, which are entry 62's own directories.
- `_p_acc_spec` now sits at `test_geo38_f_clinical_box_loaded.py:510` (item 833 cited `:376`). Import it by name.
- Count pytest verdicts per rank, by token.

**Scope.** The reader, the mapping, (K) and (D), with the degree-2 column blank. **Not claimed:**
- any degree-2 value (C37's, after entry 62);
- the move per row or the order rule (C37);
- the rod rows (C40);
- an adjudication or any AED number;
- a SPEC edit.

**Negative result.**
- (R) red means the reader leaked into the render. Fix the render's default path, keep the first log, and re-run W1 once (rule (i)).
- (K) red on the committed files means the SPEC and the page disagree as strings: a SPEC transcription question for the weekly. Report the first mismatch, edit neither file, and mark BLOCKED.
- (D-v1) red means the page's arithmetic is not the module's. That is a finding: park and BLOCKED.
- (D-v2) red means the SPEC's v2 digits do not follow from C18's constant either. Print every figure at both C values, stop, and journal it for the weekly.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 C36 → DONE. C37 is then one slot after entry 62.
- `ANS-11` stays 🟡; Done-when (1) closes on C37.
- The node-C digit moves are named for the weekly. They read as a SPEC v2 transcription question: the node is 185.7629 pF, but its C_p / Z_RC / Z_th digits were computed at 185.762904 pF, and its R_p uses rQ², not the binding ruling's r(1 + Q²). The weekly rules on them; no slot edits the SPEC.

**Depends on:** nothing. Item 880 (C43) edits the same two files and runs after this item.

**Outcome (06:00 slot, 2026-11-08): DONE, every gate green.** The page gains `read_c19`, `degree2_cells`, `comparison_cells`, `check_k`, `node_figures`, `parse_spec_node` and `check_d`; `render_comparison(m, c19=None)` takes its degree-2 cells from `degree2_cells`. W1 `20261108T120606Z_ANS-11-C36-W1.log`: 39 passed, 0 failed, 23 s (item 845's 30 plus nine new nodes). Page `20261108T120640Z_ANS-11-C36-page.log`: 64 `[PASS]`, 0 `[FAIL]`, 2 s, and the `read_c19` line reads both paths `absent`. `git diff --stat` on `metrics.json` / `COMPARISON.md` is empty, so (R) holds. Docrefs `20261108T120659Z_ANS-11-C36-docrefs.log`: `dead=0 guide=0`, exit 2 on the stale report, as item 845's. (K) is 18/18 green, with the six "not carried" tokens printed. (C19-stub) fills 16 cells, and only the copper / ESR row stays blank. (D-node) holds: 185.7629 at `:14` = `:75` = `f"{ROOT_D2_PF:.4f}"`. (D-v1) is 9/9: both `.9e` node tokens equal `:12063`, and R_p 825.0 / Q_C 64.2 / C_p 193.5527 / Z_RC / Z_th equal v1's strings. (D-v2) is 6/6 at C = `ROOT_D2_PF`. Controls (N-K), (N-D), (N-stub-key), (N-stub-node) and (N-stub-deg) are all as predicted.
- **Printed, not asserted (rule (e)), at the node C = 185.7629 pF.** R_p 896.1 and Q_C 66.9 are unchanged. C_p is 185.7214 against SPEC v2's 185.7215. The SPEC v2 strings move, in units of their last printed digit, by −Im Z_RC +289, Re Z_th +133 and −Im Z_th +23. The exact r(1 + Q²) gives 896.3 at either C. All of this is as the brief predicted, and it is named for the weekly (the SPEC v2 transcription question above). No SPEC edit.
- **Disclosed.**
  - The two docstring citations were corrected to `SPEC.md:120–126` / `:132–134`.
  - The page window ran as the brief spells it, and the page printed `FEM_EM_ANS_SKIP_PRIVATE=unset`, as item 833's page log did. `run_examples.sh` does not forward the variable, and the page has no private read path.
  - The degree-2 sense cell uses `.6f`, the C6 log's own format, so that it compares like-for-like with the degree-1 cell.

### §9 item 880 — 2026-11-08 (§10 C43, chain C — the ans:11 page gains SPEC v2's node section; lane 1, docs, -n 2)

*Written by the 2026-11-08 03:00 daily review from the 11-07 weekly's enumerated C43 (`PROJECT_PLAN.md` §10 chain C, "C43 (now; lane 1, docs, `-n 2`, no solve; serial on nothing)"), under §9 rule (4). Family C43, attempt 1 of 4. **Queued after item 871**, which adds `node_figures` / `parse_spec_node` and edits the same script and test. The weekly's "serial on nothing" holds for its inputs; the order exists only because both items touch the same files.*

**The weekly's text (binds).** "The `ans:11` page gains SPEC v2's node section (185.7629 pF; R_p = 896.1 Ω ‖ C, Q_C = 66.9, C_p = 185.7215 pF — the figures derived in the page's own code from the module's ω, C and r, asserted equal to the SPEC's by string); `docrefs` 0 dead."

**Code.** Three files only: the page script, its guide `.md`, and `tests/unit/test_ans11_record_page.py`. No SPEC edit, no `metrics.json` / `COMPARISON.md` change; (R) stays byte-identical.
1. **Script.** `node_section_lines(spec_node)` returns the section's printed lines, built only from item 871's `node_figures` / `parse_spec_node`. `main()` prints them after the `[PASS]` block:
   - the node line, which gives C = 185.7629 pF = `f"{ROOT_D2_PF:.4f}"` and its source (entry 56's crossing 185.762904 pF, `SPEC.md:17–19`);
   - one `[PASS]` / `[FAIL]` line per v2 figure (R_p, Q_C, C_p, Z_RC, Z_th) against the SPEC's string, derived at `ROOT_D2_PF`, as (D-v2);
   - one printed line of the same figures at the node C, each with its last-digit move against the SPEC's string;
   - the line "R_p by r(1 + Q²) printed beside the SPEC's high-Q form".
   A `[FAIL]` adds to `failed`, so the page exits 1.
2. **Guide.** It gains `## 4. SPEC v2's node (the capacitor the operator enters)`, after section 3. The section:
   - says the node is 185.7629 pF on all 64 gaps with r = 0.2 Ω, and that the R_p ‖ C below is HFSS's parallel entry at 64 MHz;
   - quotes **only the closing window's printed digits** (rule: closing-log digits only);
   - names the node-C digit moves as printed, not asserted, and as the weekly's to rule;
   - says the degree-1 record column above ran at 193.5996 pF, and that 185.7629 pF is the degree-2 adjudication column's node (entry 62).
   The guide's "Not claimed" paragraph drops "SPEC v2's node" and keeps "the degree-2 rows". No other guide line changes.
3. **Test.** `test_c43_node_section`:
   - every `[PASS]` line of `node_section_lines` is a PASS;
   - every numeric token in the guide's section 4 (regex `[-+−]?\d+\.\d+(?:e[-+]\d+)?`, the hyphen first so it is literal) occurs in `node_section_lines`' output, so the guide quotes no digit the page does not print;
   - the section heading occurs exactly once.

**Anchor (asserted).**
- (D-v2) through the page's printed lines: R_p 896.1, Q_C 66.9, C_p 185.7215, Z_RC −13.386935347j and Z_th 3.505604063 − 12.39875443j equal `SPEC.md:75`, `:77` by string. They are derived from `FREQUENCY_HZ_C4`, `ESR_OHM` and `ROOT_D2_PF` through `_z_rc` / `_thevenin_feed_impedance`, never retyped.
- (D-node): the SPEC's node string = `f"{ROOT_D2_PF:.4f}"`.
- The guide-digit check above.
- `docrefs` `dead=0`.

**Controls (asserted).**
- **(N-guide)** A tmp copy of the guide with one digit of the section's C_p changed fails the guide-digit check.
- **(N-node)** `node_section_lines` on a tmp SPEC whose node reads 193.5996 prints `[FAIL]` on R_p (825.0 ≠ 896.1).

**Windows.** Bare, foreground, Bash timeout 660000 ms. `pgrep -c python3` reads 0 first.
1. **`ANS-11-C43-page-pre`:** `scripts/testing/run_and_log.sh ANS-11-C43-page-pre "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"`, after the script edit and before the guide edit. It gives the printed lines the guide quotes. It closes nothing.
2. Edit the guide.
3. **`ANS-11-C43-W1`:** item 871's W1 command, labelled `ANS-11-C43-W1`, with raw log `ans11-c43-w1-raw.log`, **under `mpiexec -n 2`** (the weekly's width; every unit node is rank-agnostic or rank-0-writes-then-bcast), `timeout -k 30 180`. Expect item 871's count plus one.
4. **Closing, `ANS-11-C43-page`:** `scripts/testing/run_and_log.sh ANS-11-C43-page "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"`. This is the guide's own flags (`-n 1`; examples close on their guide's run line), run after the guide edit. Its node lines must equal the pre-window's by string; the arithmetic is deterministic. If they differ, re-quote and re-run the closing window.
5. **`ANS-11-C43-docrefs`**, item 871's form, which must read `dead=0`.
6. **`ANS-11-C43-figures`:** `scripts/testing/run_and_log.sh ANS-11-C43-figures "timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py"`. `ans:11`'s figure is `mesh:17`'s committed copy (the 11-07 weekly's ratification). Name the result.

**Tier.** Smoke / standard, `-n 2` for W1, `-n 1` for the page, complex, no solve. **Class:** ops / docs. **Cost:** ≈ 5.1 measured / 25 predicted min. Windows are measured at ≈ 3 s per page run (`20261105T223815Z_ANS-11-page.log:78`) and 24 s for W1 at `-n 1` (`20261106T013115Z_ANS-11-spec-W1.log:909`).

**Traps already paid for.**
- Item 871's list.
- A guide is never edited after its closing window. If the guide changes, re-run the closing window (the 18:00 ruling of 11-07).
- A guide quotes a round-off term as a bound, or names the window its digits came from (03:00 ruling, `EX-171`).
- Use `./scripts/run_examples.sh`, not the guide's root form.
- The census is unchanged (ans 9); this item adds no example.

**Scope.** The page's node section and the guide's text. **Not claimed:**
- any degree-2 value;
- the SPEC's node-C digits (the weekly's ruling, item 871);
- the HFSS session;
- the setup figure.

**Negative result.**
- A `[FAIL]` on a v2 figure means the SPEC's v2 strings do not follow from the module constants. Stop, BLOCKED, and journal it for the weekly. Never edit the SPEC.
- A guide-digit miss means the guide quotes a digit the page does not print. Fix the guide and re-run the closing window, once.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 C43 → DONE.
- The `ans:11` page shows the operator the node to enter, with its R_p ‖ C derived in code (`ANS-11` replication readiness).
- `ANS-11`'s glyph stays 🟡.

**Depends on:** item 871 (same script and test; uses its functions).

**Outcome (09:00 slot, 2026-11-08): DONE, every gate green.** The page gains `node_section_lines`, printed after the `[PASS]` block; a `[FAIL]` there joins `failed`. Page-pre `20261108T152453Z_ANS-11-C43-page-pre.log` and closing page `20261108T152642Z_ANS-11-C43-page.log` (2 s, exit 0, 71 `[PASS]`, 0 `[FAIL]`): the node lines are equal by string. The node line is `[PASS]` (185.7629 at `:14` = `:75` = `f"{ROOT_D2_PF:.4f}"`, the crossing 185.762904 pF). The six v2 lines are `[PASS]` at C = `ROOT_D2_PF`: R_p 896.1, Q_C 66.9, C_p 185.7215, −Im Z_RC 13.386935347, Re Z_th 3.505604063 and −Im Z_th 12.39875443, each +0 in its last digit. At the node C, printed only: R_p +0, Q_C +0, C_p 185.7214 (−1), −Im Z_RC +289, Re Z_th +133, −Im Z_th +23. r(1 + Q²) gives 896.3 at either C. (R) holds: `git status` shows no change to `metrics.json` / `COMPARISON.md`. W1 at `-n 2`: `20261108T152614Z_ANS-11-C43-W1-r2.log`, 42 passed per rank, 0 failed, 21 s. The guide's section 4 carries 15 numeric tokens and misses none. (N-guide) misses exactly `185.7216`, and (N-node) prints `[FAIL]` on R_p (824.9762, i.e. 825.0 at one decimal, vs 896.1). Docrefs `20261108T152647Z_ANS-11-C43-docrefs.log`: `dead=0 guide=0`, exit 2 on the stale report, as item 871's. Figures `20261108T152649Z_ANS-11-C43-figures.log`: 182 ok, 0 missing, 0 broken; `ans:11` ok on `mesh:17`'s copy.
- **Disclosed.**
  - The first W1, `20261108T152538Z_ANS-11-C43-W1.log` (41 passed, 1 failed per rank), is kept. Its one failure was the test's own (N-node) check: a substring `"825.0"` test against the `.7g` FAIL print `824.9762`. The control itself went red as predicted. The check now reads the printed value at the SPEC's one decimal. It was re-run once as `-r2`, and no bound changed.
  - W1 adds three nodes, not one: `test_c43_node_section`, plus the two controls `test_c43_n_guide` / `test_c43_n_node`.
  - The v2 figures are derived at `ROOT_D2_PF` when the SPEC's node names it. Otherwise they are derived at the node as printed, which is how (N-node) fails by arithmetic.
  - The guide's section 4 names the degree-1 node by its log line (`:12063`), not its digits, because the node section does not print them.

### §9 item 983 — 2026-11-14 (§10 C37, chain C — rule (4), the 11-07 weekly's step, after-window licence on `xxl` entry 62; chain step (no solve))

*Written by the 2026-11-14 10:30 daily review from the 11-07 weekly's C37 (`PROJECT_PLAN.md` §10 chain C: "`ANS-11`'s degree-2 rows filled from entry 62's JSONs, the degree-1 → degree-2 move per row printed, the 10-24 exit's order rule applied; JSON written. A record; `ANS-11`'s Done-when (1) closes on it."), under §9 rule (4). It does not count against rule (6). Family C37, attempt 1 of 4. **Takeable now under the after-window licence:** this review read entry 62 mechanically and wrote its `xxl-ledger.md` row (2026-11-14 10:30, from a `log-pathologist` reading). Entry 62's decision rule's green branch *is* this step ("`ANS-11`'s `COMPARISON.md` rows are filled from it, the moves are records, the order rule applies"); what the order rule's verdict then means for the clinical figures is the weekly's.*

**One correction of route, not of target (binding ruling, §9 digest: "every new step's reader of a pending window parses its tracked log, not the host JSON … the host JSON, where present, is only a custody check").** The weekly wrote "from entry 62's JSONs" before that ruling. The JSONs are gitignored (`/workspace/logs/geo38/c6-node-185.7629-d2.json`, `/workspace/logs/geo36/i36-node-185.7629-d2.json`), so a `COMPARISON.md` rendered from them cannot regenerate on a clone. **The degree-2 cells are therefore parsed from entry 62's tracked log, `docs/testing/logs/20261114T080009Z_GEO-38-C19.log`, and C36's `read_c19()` JSONs become the custody check (J62).** The rows, the moves and the order rule are the weekly's, unchanged.

**Facts verified for this brief (10:30 review, by `grep` on the tracked log).**
- Every `_SINGLETONS` line of the page prints **exactly once** in entry 62's log, all in call 1 (call 1 ends at `[C19] call 1 … rc1=0`, `:12703`; call 2's lines carry other tags): `(Z) node 185.7629 pF` `:12160`; `kinds … S2` `:12379`; `Y_in ccw` `:12383`; `Y_in cw` `:12384`; `mean |B1+|_ccw` `:12387`; `c*` `:12406`; `normalisers` `:12410`; `wrote …c6-node-185.7629-d2.json` `:12411`; `price` `:12412`. The 18 `lattice` lines are present (count 18). So **`parse_c6_log(C19_LOG)` reads entry 62 unedited** — its own duplicate check is the "exactly one distinct match" rule (10:30 11-12 ruling).
- The box-wide grid's result is one line, `:12407`: `[GEO-38 C6a] (C6b) grid n 90 (72 new balls, 18 lattice values reused); g* x+157.5mm y+0.0mm dz+0.0mm (0.157500, 0.000000, 0.000000) m: SAR10g <F> W/kg (ccw), (1,+1j) ball <F> W/kg; beside c* …`. Its label form (`x… y… dz…`, one decimal) is **not** `_LABEL`'s; it needs its own regex. `:12162` is the grid's census line (`n 90 of 91 candidates kept …`).
- The `(Z) node` line ends `… failures none (ASSERTED)`; `_SINGLETONS["node"]` uses `re.search`, so the suffix does not matter.
- **Entry 62 is green on every asserted gate** (`log-pathologist` reading, ledger row 2026-11-14): Status 0, `rc1=0` `:12703`, `rc2=0` `:25454`, 16 / 16 rank summaries per call, FAILED / ERROR 0. **The log prints no degree-1 comparand**, so every move is computed by this step from entry 52's tracked log (the page's `SOURCE_LOG`), never read from entry 62. In call 1 "P_rod" (`:12380`) is tag 403 counted as gel, not a rod. Whether g\* at x = +157.5 mm lies on the grid's edge is not established: print g\*'s coordinates beside the grid's extent if the module exposes it, and claim nothing about a maximum outside the grid.

**Code: the `ans:11` page script and `tests/unit/test_ans11_record_page.py`; the page window regenerates `ans:11`'s `metrics.json` and `COMPARISON.md`. Nothing else** (no `src/`, no SPEC, no C6 / I36 / C33 / `ans:12` / `ans:13` module edit).
1. **`C19_LOG_BASENAME = "20261114T080009Z_GEO-38-C19.log"`, `C19_LOG`** beside `SOURCE_LOG` (read through `_read_log`, so a gzipped copy still reads).
2. **`parse_c19_grid(path)`** — one new regex for `:12407`'s form; returns `{"line", "n", "label", "x_m", "y_m", "z_m", "sar10g", "sar10g_cw"}` as printed strings; raises `ValueError` on zero or more than one match.
3. **`degree2_block(path=C19_LOG)`** — `p2 = parse_c6_log(path)`, `d2 = derived(p2)`, `g = parse_c19_grid(path)`; asserts `p2["node"]["c_pf"] == C19_NODE_PF`; returns `{"source_log": C19_LOG_BASENAME, "parsed": p2, "derived": <the same derived-key dict `build_metrics` writes, from d2>, "grid": g, "moves": …, "order_rule": …}`.
   - **`moves`**: for each row with a degree-1 and a degree-2 number — S₁₁ / S₂₂ magnitude, |Y_in| ccw P33 / P41, the lattice mean and its per-√P_acc value, SAR10g(c\*) and per P_acc, the cw ball, the centroid mean and its per-√P_acc value, P_box / P_src / P_acc, the gel fraction, the sense — `f"{x2 / x1 - 1:+.3e}"` (≤ 4 significant digits; a ratio, never an AED figure).
   - **`order_rule`**: the 10-24 exit's rule applied to the **per-watt rows** (lattice mean per √P_acc, SAR10g(c\*) per P_acc, centroid mean per √P_acc, gel fraction P_box / P_acc): the list of rows whose |move| > 2 %, and the sentence "clinical figures quoted at degree 2 only" if that list is non-empty, else "the 2 % order rule does not fire". **A printed verdict, not an assertion**; the 2 % is the exit's pre-registered threshold, not a band this step sets.
   - c\*'s label moving (degree-1 c\* vs degree-2 c\*) is printed as a fact; the box-wide g\* is a degree-2-only row (entry 52 ran no grid).
4. **`build_metrics(parsed, d, case="ans11")`**: for `"ans11"`, `"degree2_c19"` becomes `degree2_block()` (the tracked log, so `build_metrics(parsed, derived(parsed))` stays a pure function of tracked files and the landed nodes `test_r_regeneration` / `test_c49_k_ans11_default_byte_identical` stay green **unedited** on the regenerated files). `"ans12"` is unchanged (`"pending (entry 66)"`); ANS-12's `COMPARISON.md` byte-identical.
5. **`render_comparison(m, c19=None, case)`**: when `c19` is given (C36's stub path) the cells come from `degree2_cells(c19)` as today, so `test_c19_stub` stays green unedited; **when `c19` is `None` and `m["degree2_c19"]` is a dict, the degree-2 cells come from that block**, in the **same cell formats as the degree-1 column** (row by row, the d1 f-string with the degree-2 tokens and their `:line`), so each row reads like-for-like. The grid row gets `g`'s label and SAR10g with `:12407`. The copper / ESR row stays `""` (entry 62 prints no per-partition fractions either; say so in its cell text only if `CASES` already carries such text for d1). A third section **"Degree 1 → degree 2 (records)"** lists `moves` and the `order_rule` sentence. `CASES["ans11"]`'s `title`, `adjudication`, `h2`, `not_claimed` texts are updated to say the degree-2 column is entry 62's record (`xxl` 62, `-n 16`, degree 2) and still un-adjudicated (Done-when (2) is the weekly's). **Before editing `CASES["ans11"]`, `grep -rn 'CASES\[' tests/ examples/`**: any reader found must stay green unedited, or stop.
6. **`main()`** passes `None` (not `read_c19()`) to `render_comparison`, so the committed render is the tracked-log render; it still calls `read_c19()` and runs **(J62)** below, printing one `[PASS]`/`[FAIL]` line; and it prints the `moves` table and the `order_rule` sentence.
7. **`custody_c19(block, c19)`** — (J62): for each raw field the C6 JSON carries and the log prints (S2 four entries, Y_in ccw / cw P33 / P41, the 18 lattice `b_plus_ccw` / SAR values by label, c\*'s label / SAR / cw SAR, `mean_b_plus["ccw"]`, P_box / P_src / P_acc, the grid's star label / SAR), format the JSON value **in the log token's own print format** (`.10e` for S, `.9e` elsewhere, the label as printed) and compare strings. Returns `(n_fields, mismatches, print_form_only)`. Per the 10:30 11-12 ruling a value-equal, print-form-only difference is listed in `print_form_only` and printed, **not** red; a value mismatch is red. With the JSON absent (a clone) it returns `None` and prints "absent; custody not run" — the tracked-log route is asserted regardless. The I36 JSON is read and validated by `read_c19` only (no row; C40's family reads the rod).

**Gates (asserted, in the unit test and the page's `[PASS]` lines).**
- **(X2)** every degree-2 cell of the regenerated `COMPARISON.md` contains its row's entry-62 token string (via `comparison_cells`' degree-2 twin), with its `:line`; the lattice table's 18 degree-2 rows, if rendered, likewise.
- **(D2)** the log-internal identities, each at **the tolerance `check_derived` already uses for the same quantity** (imported, not new): `_p_acc_spec(S, (1, −1j), V_s, 50)` from `:12379` / `:12160` vs `:12410`'s printed P_acc (rel 1e-8); SAR10g(c\*) / P_acc vs `:12410`'s `sar_per_p_acc` token (5e-10); P_box / P_acc vs `:12410`'s printed ratio (5e-9); the centroid mean / √P_acc vs `:12410`'s `b1_per_sqrt_p_acc` token at its printed digits; argmax(lattice SAR10g) == c\*'s label; c\*'s SAR string == its lattice line's string. **Not** `check_derived` itself (its SPEC comparisons are degree-1 tokens). **Trap, read at this review: entry 62's S₁₂ and S₂₁ strings differ in the 11th significant digit** (`:12379`: `4.2062703482e-03+9.9390915599e-02j` vs `4.2062703483e-03+9.9390915598e-02j`), so `check_derived`'s "(D) S12 string == S21 string" would be red at degree 2. Print both, never assert string equality; reciprocity is (R1), asserted in-window at 2.415e-13 (`:12379`). The degree-2 lattice mean is arithmetic on `:12388–12405`, as at degree 1 (the log prints none).
- **(J62)** with both JSONs present on this host: `mismatches == []`; `n_fields` printed (expect 4 + 4 + 36 + 3 + 1 + 3 + 2 = 53; print the actual and go on if the brief's arithmetic is off by a field the JSON does not carry — the gate is zero value mismatches, not the count).
- **(R)** `metrics.json` and `COMPARISON.md` regenerate byte for byte from tracked files (landed `test_r_regeneration`, `test_c49_k_ans11_default_byte_identical`, unedited); after the closing windows `git diff --stat` on both is empty and ANS-12's `COMPARISON.md` is unchanged.
- **(K)** C36's `check_k` green on the new `COMPARISON.md` (degree-1 column unchanged, the six "not carried" tokens printed).
- **Every existing node green:** `test_ans11_record_page.py`, `test_ans12_record_page.py`, `test_ans13_record_page.py`, `test_doc_reference_exit_codes.py`; C33's module (knob `1`) and P89′'s module as item 979 ran them; `-e ans:12`, `-e ans:13` ((R-cmp)), `-e mri:41`.
- **`test_c19_absent` is retired** — its premise ("entry 62 has not run") expired at `3185acfd` (known-issues 2026-11-14) — and **replaced by `test_c37_custody`**: (J62) as above when present, and the log route asserted always. Disclose the retirement in the commit (a premise node, not a gate; the item-918 precedent).
- `docrefs` `dead=0`.

**Controls (asserted).**
- **(N-J62)** a tmp copy of the C6 JSON with one lattice `b_plus_ccw` value scaled by `1 + 1e-6` makes (J62) red on that label alone (the `.9e` print resolves 1e-9 relative, so the ceiling is ample).
- **(N-X2)** a tmp copy of entry 62's log with one digit of `:12410`'s P_acc changed makes (D2) red (and (X2) on the P row).
- **(N-grid)** a tmp log with `:12407` duplicated makes `parse_c19_grid` raise.
- **(N-node)** a tmp log whose `(Z) node` reads `193.5996` makes `degree2_block` raise.
- C36's / C43's / C49's controls unchanged.

**Printed (records):** the 17 degree-2 rows, the moves table, the order-rule sentence, (J62)'s `n_fields` and any `print_form_only` entries, g\* beside c\*.

**Windows.** Bare, foreground, Bash `timeout: 660000`; `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 first. `FEM_EM_ANS_SKIP_PRIVATE=1` on every window.
0. **Pre-check:** `xxl-ledger.md`'s 2026-11-14 `GEO-38-C19` row is filled (this review); `git status` clean; no other slot this interval has edited the page.
1. Code.
2. **Page (writes the files):** `scripts/testing/run_and_log.sh ANS-11-C37-page "./scripts/run_examples.sh -e ans:11 -n 1 -t 180"` — every `[PASS]`, (J62) `[PASS]`, the moves and order-rule lines printed. `git diff --stat` now shows exactly `metrics.json` and `COMPARISON.md` of `ans:11` changed (and nothing under `ans12_…`).
3. **W1:** C49's W1 command (`20261110T104946Z_ANS-12-C49-W1.log:12`, `mpiexec -n 2`, `timeout -k 30 180`) with label `ANS-11-C37-W1`, raw log `ans11-c37-w1-raw.log`, paths `tests/environment/test_complex_mode.py tests/unit/test_ans11_record_page.py tests/unit/test_ans12_record_page.py tests/unit/test_ans13_record_page.py tests/unit/test_doc_reference_exit_codes.py`. Expect 0 failed per rank (count by token).
4. **C33 / P89′ regressions:** item 979's two commands verbatim (`20261114T133421Z_GEO-38-C33c-W1.log:12`, `20261114T133435Z_MAT-21-P89-W1c.log:12`), labels `ANS-11-C37-c33` / `ANS-11-C37-p89`.
5. `scripts/testing/run_and_log.sh ANS-11-C37-ans12-page "./scripts/run_examples.sh -e ans:12 -n 1 -t 180"`, then `git diff --stat` on ANS-12's `COMPARISON.md` empty; `scripts/testing/run_and_log.sh ANS-11-C37-ans13-page "./scripts/run_examples.sh -e ans:13 -n 1 -t 180"`; `scripts/testing/run_and_log.sh ANS-11-C37-mri41 "./scripts/run_examples.sh -e mri:41 -n 2 -t 120"`.
6. **Closing page re-run** (same command, label `ANS-11-C37-page-r`): `git diff --stat` on `ans:11`'s two files is now empty against the files written in window 2 (R). This is the closing window of Done-when (1).
7. **Docrefs:** item 871's form, label `ANS-11-C37-docrefs`, `dead=0`.

**Tier.** Standard, `-n 1` (W1 and the regressions at `-n 2`), complex (the C6 import), no solve. **Class:** chain step (no solve). **Cost:** ≈ 10.8 measured / 40 predicted min; the windows sum to ≈ 3 min (C36's 23 s W1, C49's 2 s pages, item 979's two windows ≈ 30 s).

**Traps already paid for.**
- Item 871's / 931's / 959's lists: `re.search`, never `startswith`; either log suffix; derived values at ≤ 7 significant digits (the leak hook's floor is 8 — raw tokens are the tracked log's own and already passed it); `FEM_EM_ANS_SKIP_PRIVATE=1`; open no `aed_results/`, `COMPARISON_private.md` or `docs/private/` file; never `--no-verify`; on a hook refusal move the log unchanged to `logs/leak-held/` and report.
- **No edit to `parse_c6_log`, `derived`, `parse_spec`, `check_k` or `CASES["ans12"]`:** C33, P89′, `mri:41` and the `ans:12` / `ans:13` pages read them.
- **No top-level import added** (the page's import closure is C33's and `mri:41`'s); the grid regex is local.
- **The default render must not read a gitignored file** — that is the defect item 959 found (`main` passed `read_c19()` to the render). (J62) reads the JSON; the render never does.
- **Never write under `/workspace/logs/geo38` or `/workspace/logs/geo36`** — entry 62's JSONs are the custody evidence. Stubs and tmp copies live on `tmp_path` (rank 0 writes, then a bcast).
- **Rule (i):** no page or test edit after its last window; if the page changes after window 2, re-run 2–6.
- A changed shared constant is an edit to every module that asserts it (the 03:00 11-14 trap): grep `CASES[`, `METRICS_JSON`, `COMPARISON_MD`, `render_comparison(` and `build_metrics(` across `tests/` and `examples/` before editing, and list the hits in the outcome.
- Count pytest verdicts per rank, by token. Never delete a log, red or not.

**Scope.** `ANS-11`'s degree-2 record column from entry 62's tracked log, its moves and the order rule's printed verdict, with JSON custody. **Not claimed:** any adjudication or AED number (Done-when (2), the weekly's, with the operator's First Order column); the rod rows (entry 62's call 2; C40's family); any change to the SPEC or its node; what the order rule's verdict means for other clinical figures (the weekly's); the 3 T twin (entry 66).

**Negative result.**
- (J62) value mismatch ⇒ the tracked log and the host JSON disagree for one window — a custody finding: stop, BLOCKED, known-issues; render nothing new.
- (D2) red ⇒ the page's P_acc arithmetic does not reproduce the module's degree-2 print — a finding (the identity held at degree 1): park, BLOCKED.
- Any landed node red (C33, P89′, `ans:12`, `ans:13`, `mri:41`, the ANS-12 file moving) ⇒ stop, BLOCKED; known-issues is the review's call.
- `parse_c6_log` raising on entry 62 ⇒ the "exactly once" fact above is wrong; stop and report the line.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 C37 → DONE; **`ANS-11` Done-when (1) met** (an `ans:11` page that imports C6's module, `COMPARISON.md` with our degree-1 and degree-2 rows and the AED columns blank, closing window `./run_examples.sh`). The glyph stays 🟡 (Done-when (2) is the adjudication); the §7 state line records (1).
- Known-issues 2026-11-14 retired (`test_c19_absent` replaced; the default render no longer reads the host JSON), which **unblocks item 959's re-pose (item 984)**.
- The moves and the order-rule verdict are records for the 11-16 weekly (entry 62's decision rule names it).

**Depends on:** nothing open (entry 62 has run and its row is written). **Never in the same slot as item 984** (both edit the page).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 983 — `ANS-11` C37 — per `docs/planning/chunks/ANS-11.md` "§9 item 983". First confirm the 2026-11-14 `GEO-38-C19` row in `docs/testing/xxl-ledger.md` is filled; otherwise stop. Edit only the `ans:11` page script and `tests/unit/test_ans11_record_page.py`; `metrics.json` / `COMPARISON.md` of `ans:11` change only through the page window. Parse entry 62's tracked log for the degree-2 column; the host JSONs are custody (J62) only, and the default render never reads them. Do not edit `parse_c6_log`, `derived`, `parse_spec`, `check_k` or `CASES["ans12"]`; ANS-12's `COMPARISON.md` stays byte-identical. Retire `test_c19_absent` in favour of `test_c37_custody`, disclosed. `FEM_EM_ANS_SKIP_PRIVATE=1`; open nothing under `aed_results/` or `docs/private/`. No solve. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (12:00 CST slot, 2026-11-14): DONE, every gate green; `ANS-11` Done-when (1) met.** The pre-check passed: the ledger row is filled, the tree was clean, and no other slot touched the page this interval. The page gains `C19_LOG`, `parse_c19_grid` (its own regex for `:12407`), `degree2_block`, `moves_and_rule`, `check_d2`, `x2_misses` and `custody_c19`. Two refactors are pure: `_derived_strings` (out of `build_metrics`) and `_column_cells` (out of `render_comparison`), and the degree-1 cells stay byte-identical. `parse_c6_log`, `derived`, `parse_spec`, `check_k` and `CASES["ans12"]` are unedited. `build_metrics(…, "ans11")` now sets `"degree2_c19"` to `degree2_block()`. `render_comparison(m, None)` fills the degree-2 cells from that block, adds a 7-column degree-2 lattice table (7 cells, so `comparison_cells` ignores it and (K) is unchanged) and a "Degree 1 → degree 2 (records)" section, and names the degree-2 source in Provenance. `main()` renders with `None`. It still runs `read_c19()`, but only for (J62). The `CASES[` / `METRICS_JSON` / `COMPARISON_MD` / `render_comparison(` / `build_metrics(` hits outside the page are all in `test_ans11_record_page.py` and the `ans:12` / `ans:13` pages, and every one is green unedited. **`test_c19_absent` is retired** (its premise expired at `3185acfd`) and replaced by `test_c37_custody`. Other new nodes: `test_c37_lines`, `_d2`, `_x2`, `_moves_rule`, and the controls `_n_j62`, `_n_x2`, `_n_grid`, `_n_node`.

The windows, in order:
- Page `20261114T181052Z_ANS-11-C37-page.log`: 80 `[PASS]`, 0 `[FAIL]`, Status 0.
- W1 `20261114T181121Z_ANS-11-C37-W1.log` (`-n 2`): 76 passed per rank, 0 failed / skipped, `[capture] rc=0`.
- C33 `…181201Z_ANS-11-C37-c33.log`: 15 passed per rank. P89′ `…181206Z_ANS-11-C37-p89.log`: 9 passed. Both match item 979's counts.
- `ans:12` `…181222Z_ANS-11-C37-ans12-page.log`: (R-cmp) 5051 B == 5051 B, and `git diff` on ANS-12's files is empty.
- `ans:13` `…181228Z_ANS-11-C37-ans13-page.log`: Status 0. `mri:41` `…181231Z_ANS-11-C37-mri41.log`: Status 0.
- **Closing** `20261114T181245Z_ANS-11-C37-page-r.log`: `git diff` on `ans:11`'s two files is empty against window 2's, so (R) holds.
- Docrefs `…181259Z_ANS-11-C37-docrefs.log`: `dead=0 guide=0`, exit 2 on the stale report, as in item 871.

The page windows ran the brief's literal command, which leaves `FEM_EM_ANS_SKIP_PRIVATE` unset. The page prints `=unset`, as C36 / C49 did. W1 set it to `1`. The page has no private side either way.

**(J62): n_fields 53, mismatches 0, print_form_only none.** **(D2):** the `_p_acc_spec` identity holds at rel 1.751e-10, SAR/P_acc at 6.9e-11, P_box/P_acc at 1.1e-10, and centroid mean / √P_acc at 1.0e-15 against a half-ulp of 5e-15. The argmax is c\*, and the c\* string equals its lattice line's string. The degree-2 S₁₂ and S₂₁ strings differ in the last digit; that is printed, not asserted. **(X2)** misses: none.

**Moves (degree 2 / degree 1 − 1):** |S₁₁| −1.129e-01, |S₂₂| −9.305e-02, |Y_in| P33 −1.121e-01, P41 −1.375e-01, lattice mean −1.272e-02, **per √P_acc −4.624e-03**, SAR10g(c\*) −2.349e-02, **per P_acc −7.412e-03**, cw ball −5.164e-02, centroid mean −1.302e-02, **per √P_acc −4.924e-03**, P_box −2.859e-02, P_src +1.772e-02, P_acc −1.620e-02, **gel fraction −1.259e-02**, sense −3.4e-07. **Order rule: no per-watt row exceeds 2 % → "the 2 % order rule does not fire"** (a printed verdict). c\* stays x−105mm dz−50mm. The box-wide g\* is x+157.5mm y+0.0mm dz+0.0mm. The host JSON's grid spans x [−157.5, +157.5] mm, so **g\* sits on the grid's outermost x plane** (custody print; no claim about a maximum outside the grid). **Caveat, not in the brief:** entry 52 ran at C = 193.5996 pF and entry 62 at 185.7629 pF, so every move carries the node change as well as the order change. `COMPARISON.md`'s not-claimed line says so. Whether that bears on the order rule's reading is the 11-16 weekly's call.
