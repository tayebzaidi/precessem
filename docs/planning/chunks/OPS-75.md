# OPS-75 — the two example scripts that still read a dated log by its plain name

*History file for the §7 `OPS-75` row (`OPS-46`). Opened by the 2026-11-01 18:00 daily review, lane 1, from `OPS-73`'s disclosure (3): "`H_TABLE_LOGS` in `04_birdcage_four_port_10_64_128MHz.py` and `node_log_readings` in `08_f_human_quadrature_mode_content.py` still read plain names only, and the tests go round them with scratch copies, so any route that runs them directly on those logs fails."*

**Done-when.** Both scripts read a record whether retention left it plain or `.gz`; (H′) and (A′) below asserted; the tests' scratch copies stay, untouched.

### §9 item 751 — 2026-11-01 (lane 1, `OPS-73`'s disclosure)

**Why.** The same defect `OPS-73` fixed in `tests/`, on the one route it left: running the example. An example cannot import `tests/`, so the fallback is written in the script, stdlib only.

**Anchor.** New `tests/unit/test_example_log_reads.py`. **(H′)** `node_log_readings` on a `.gz`-only copy of a log returns the same dict as on the plain copy, and `H_TABLE_LOGS`' reader the same text, equal to `tests/log_records.read_log_text` on that path. **(A′)** with neither file present each raises `FileNotFoundError`.

**Negative control.** The script's read as it stands on `main` raises on a `.gz`-only path: pin this as `open(path)` on that file raising, in the test, as `test_log_records.py:69–72` does.

**Execution.** A small `_read_log(path)` in each script, or one in each script's shared helper if there is one: `gzip.open(…, "rt", errors="replace")` when `path` is absent and `path + ".gz"` exists. `08_…py:426` (`with open(log_path) as fh:`) and `04_…py`'s reads named by `H_TABLE_LOGS` (`:289`). Load the scripts as `tests/unit/test_ans8_v3_d2_block.py:_load` does. Build fixtures under `tmp_path` from short synthetic text.

**Tier, ranks, wall clock.** Standard, `-n 1`, `timeout -k 30 170`; ≈ 5 s. Then the two modules whose scratch copies touch these reads, `tests/unit/test_ans4_h_table.py` and `tests/validation/test_ans8_v3_record_rows.py`, as their `OPS-73` windows ran. **Class:** ops / docs. **Cost:** ≈ 5.6 measured / 15 predicted min. **Depends on:** nothing; item 740 edits the same `08_…py`, so if a slot takes both, take 740 first.

**Traps already paid for.** **Never `gunzip` a record to get green.** A resolved `.gz` name broke (X) in `test_ans4_h_table.py` because the numbers file records each log's listed name (`20261101T182235Z_OPS-73-Wreg-h.log`): the return value must not carry the `.gz` name. Both scripts run in the complex container, and no change may touch a solve path. `-n 1` for a module that is not rank-safe by design.

**Scope.** Two read routes. **Not claimed:** retention, a record, a value. The example scripts' other reads are JSON, not logs, and stay.

**Negative result.** A log that is gone in both forms stops that read: name it, BLOCKED if it is on the example's main path. No `PREDICTED RED:` line.

**Status it can move.** `OPS-73`'s disclosure (3) closes; `OPS-75` ✅.

**Outcome (2026-11-02, 00:00 slot): complete, `OPS-75` ✅, run inline.** A `_read_log` in each script: `04_…py`'s `parse_h_table_logs` and `08_…py`'s `node_log_readings` read through it, and `rec["log"]` stays `p.name`, the listed name. New `tests/unit/test_example_log_reads.py`. W1 `20261102T062832Z_OPS-75-W1.log` (`-n 1`, complex, Status 0, 3 s): 10 passed, the module as committed. (H′) on `.gz`-only synthetic copies: the node dict and the h rung equal the plain copy's, the listed name kept, and each `_read_log` equals `read_log_text`. (A′) `FileNotFoundError` from both. Control: bare `open` raises on the `.gz`-only path. Regressions: `…062812Z_OPS-75-Wreg.log` (`-n 1`: `test_ans4_h_table.py`, `test_ans8_v3_record_rows.py`, `test_ans8_v3_d2_block.py`, `test_log_records.py`, 9 passed, 3 skipped), `…062835Z_OPS-75-Wreg-n2.log` (`-n 2`: the two ANS-8 modules plus `test_complex_mode.py`, 7 passed, 3 skipped per rank), both Status 0.
- **Beyond the brief, one test:** `H_TABLE_LOGS` read directly on the real logs, no scratch copy, reproduces `h_table_s.json` byte for byte, with the listed names on disk `[True, False, False]` (`:59`). That is the disclosed route on real records.
- **Disclosed:** plain reads now take `errors="replace"`, as `read_log_text` does, and a `.gz` `path` still goes through `gzip`, so `H_TABLE_LOGS`' ×1 rung is unchanged (tested). The 3 skips are the knob-gated record nodes of `test_ans8_v3_record_rows.py`, which read gitignored JSON, so `node_log_readings` has run on synthetic logs here, not on a real node log. `…062804Z_OPS-75-W1.log` is the module before the extra test, kept. `Wreg` predates that test-only edit; neither script changed after it. The tests' scratch copies stay, untouched. No solve path touched, neither example run.
