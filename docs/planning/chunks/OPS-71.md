# OPS-71 — re-record `mri:2`'s guide table to the 0.11 image

Opened 2026-10-19 18:00 daily review (§9 rule (5) lane 1, a stale record found by an executed window). Item 549's journal (`POST-8` 3b-ii-e1b, 12:00 slot) noted that `mri:2`'s re-gate `20261019T170830Z_POST-8-step3b-ii-e1b-mri2-regate.log` (Status 0, 10 s) prints digits its guide does not carry. The guide's "On record" table is the 2026-08-08 run (`20260808T020414Z_EX-3-gate.log`), from before the 0.11 image. `mri:2` does not use the sampler, so the drift is not `POST-8`'s. It is mesher-dependent record drift, which is the `(1*)` licence (`OPS-31` / `OPS-33` precedent): re-record the record, move no band.

Measured site list (guide `examples/mri/02_mass_averaged_sar.md` → the 17:08 log):

| Guide site | Guide reads | Log reads |
|---|---|---|
| table `:56`, mesh cells | 74 216 | 74 218 (`:114`) |
| table `:50` and block `:94`, DG0 sphere-average vs closed form | 1.32e-15 | 2.81e-15 (`:119`) |
| table `:54` and block `:120`, 10 g kernel mass | 0.0044 % | 0.0041 % (`:128`); also 0.0041 % in `20261018T125346Z_POST-8-step3b-W5a-regate.log:223` |
| table `:56` / `:82–83`, mesh time | 7.1 s, 14 s wall | 2.7 s, 10 s wall (timings, informational) |

Every other row matches: point 4.96e-16 (`:118`), 1 g ratio and kernel 0.0120 % (`:123–124`), 10 g ratio (`:127`), separation 2.1894 against 1/f = 2.1681 (`:131`).

### §9 item 552 — 2026-10-19 (lane 1; ops / docs)

**Executor:** `record-reconciler`, spawned foreground. **Tier:** standard, `-n 2`, one example window of ≈ 10 s plus two census windows. **Cost:** ≈ 3 measured minutes (ops / docs class) against 15 predicted. **Depends on:** nothing. Any slot may take it.

*No review ruling is open on this item.*

1. **Anchor.** A fresh run of record of the **unedited** `examples/mri/02_mass_averaged_sar.py` at `-n 2`, complex, must pass all of the example's own assertions at their unchanged bands:
   - closed-form round-off ≤ 1e-12 on both paths;
   - the averaging identity within 0.5 % at 1 g and 10 g;
   - kernel mass within 0.1 %;
   - the separation control floor 1.5 and ceiling ±5 %.

   Then every digit in the guide's "On record" table (`:42–56`) and its printed-output blocks (`:91–121`) must equal that log's printed digit. The journal holds a line-by-line trace: guide line → log line.
2. **Negative control (asserted, rule (e)).** Before the edit, three rows already differ from the 17:08 log on the same script and image: cells, DG0 and the 10 g kernel (table above, with log lines). The executor first confirms that the pre-edit guide still carries 74 216 / 1.32e-15 / 0.0044 %. After the edit, the trace shows zero mismatches. A trace that finds no pre-edit mismatch means the guide was already fixed: stop and journal.
3. **Commands** (each is the whole Bash command, bare and repo-relative; `timeout: 660000`):
   - `scripts/testing/run_and_log.sh OPS-71-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`. This is the `Command:` line of `20261019T171453Z_EX-136-pre-docrefs.log:12`.
   - `scripts/testing/run_and_log.sh OPS-71-mri2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 examples/mri/02_mass_averaged_sar.py'"`. This is the 17:08 log's `Command:` line verbatim.
   - The post-edit docrefs census, the same command as the pre-edit one, with the chunk tag `OPS-71-post-docrefs`.
4. **Traps already paid for.**
   - The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. **The spawning slot puts that sentence in the spawn prompt verbatim.**
   - Never edit `02_mass_averaged_sar.py`; it is the gate's twin (`MAT-4` step 3).
   - Never touch a band, a quadrature degree or a budget sentence (`:58–59`).
   - Docrefs gates on `exit != 1`; exit 2 (stale) is information.
   - The `(1*)` licence covers records only.
5. **Scope boundary.** `mri:2`'s guide only. This item does not sweep other example guides for the same drift, does not touch `MAT-4`'s gate records, and is not a gate change. Timings (`:56`, `:82–83`) are updated to the run of record as information.
6. **Negative result.**
   - **Any example assertion red:** stop, journal the reading, and open a known-issues entry. It would be a regression since 17:08 today on an unchanged script. Re-record nothing.
   - **The DG0 round-off reads a value other than 2.81e-15:** it is round-off and varies run to run. Write that row as "round-off, ≤ 1e-12 (2.81e-15 and <new> on <both logs>)", not as a single digit, and say so in the journal.
   - **The cell count differs from 74 218:** the mesh is not reproducible on this image. Stop, journal both counts, and hand to the review.
7. **Status it can move.** `OPS-71` ⬜ → ✅. `mri:2`'s guide records point to a run of record on the current image, with the 2026-08-08 run kept as the historical line ("was 74 216 / 1.32e-15 / 0.0044 % on `20260808T020414Z_EX-3-gate.log`").

**Commit:** the guide, the three logs, their `test-results.md` rows, the §7 row (✅, with the run-of-record log), the §9 line collapsed to DONE, and an outcome line here.

**Spawn prompt (record-reconciler):** "Execute `OPS-71`, §9 item 552, from `docs/planning/chunks/OPS-71.md` § '§9 item 552'. Measured site list: the table at the top of that file. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run every harness window in the foreground with Bash `timeout: 660000`. Never edit the example script or any band. Return the guide line → log line trace and the three log paths."

### §9 item 552

**Outcome (2026-10-20, 19:30 take-next slot, record-reconciler): DONE.** Run of record `docs/testing/logs/20261020T004018Z_OPS-71-mri2.log` (Status 0, 10 s harness, 8.3 s internal; unedited script, all assertions green at unchanged bands). Pre-edit negative control confirmed (guide carried 74 216 / 1.32e-15 / 0.0044 %). Re-recorded: cells 74 218 (log `:114`), DG0 2.81e-15 (`:119`, equals the 17:08 value so no round-off-range wording), 10 g kernel 0.0041 % (`:128`), mesh 2.6 s, 10 s wall, 8.3 s internal. Trace (guide line → log line): table 50 / block 94 → `:119`; table 54 / block 120 → `:128`; table 56 → `:114`; table 49 / block 93 → `:118` (unchanged 4.96e-16); table 53 / block 119 → `:124`; table 55 / block 131 → `:131`; closed form → `:116`; ratios → `:123`, `:127`. Docrefs census exit 2 (118 unrelated stale artifacts) before (`20261020T004014Z_OPS-71-pre-docrefs.log`) and after (`20261020T004041Z_OPS-71-post-docrefs.log`), identical but for ages; both `exit != 1`. Collected-count identity: not applicable (example, not pytest).
