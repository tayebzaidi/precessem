# OPS-80 — which callers of the two-tag power-identity helper run on a mesh with other lossy cells

*History file for the §7 `OPS-80` row (`OPS-46`). Opened by the 2026-11-04 18:00 daily review under §9 rule (5), lane 1: an open known-issues entry whose own "Resolves with" a slot can execute (known-issues "OPEN 2026-11-03 (`GEO-36` I57, §9 item 801 …): the birdcage power-identity helper sums volume loss over tags 3 and 1 only").*

**The question.** `_loss_power_w` (`tests/validation/test_port_drive_superposition.py:205–221`), which `_exact_shares` uses, integrates ½σ|E|² over `PHANTOM_CELL_TAG` (3) and `CONDUCTOR_CELL_TAG` (1) only. On I41's kept-rod build that left out tag 403's loss, and (P) read 2.882e-06 against 1e-6 (`20261103T195044Z_GEO-36-I57-diag-W1.log:1919–1942`). The entry's Scope says "the callers have not been audited for this", and its Resolves-with names that audit beside the weekly's call on widening the helper. **This item is the audit, and nothing else.** At HEAD `b29c4f77`, `grep -rnw '_exact_shares\|_loss_power_w' tests src examples scripts` returns 46 lines in 16 files. `src/fem_em_solver/ports/shares.py` mentions the name in its docstring only (`:4`).

**Done-when.** A table of every `grep -rnw` line for the two names, each classified as one of: definition, import, prose/docstring, or a call. Each call is classified on the mesh it runs on:
- **(U) unaffected:** the only lossy volume cells are tags 3 and 1. This covers hole builds, where tag 402 is a facet.
- **(K) kept lossy cells outside 3 and 1:** for example `keep_cells=True` with tag 403 mapped to a lossy material. It is marked **two-tag** if the gate's (P) uses the helper's sum, or **all-cells** if the caller sums over all cells.
- **(?) undecided by code reading.**

The count equals `grep`'s per name, with 0 lines unclassified. The known-issues entry's Scope gains the table's result. No test, helper, band or `src/` edit is made.

### §9 item 817 — 2026-11-04 (lane 1, known-issues "Resolves with"; ops / docs)

**Licence.** §9 rule (5), lane 1: the entry's own Resolves-with names this audit. No target is set. **Whether the shared helper is widened stays the weekly's call** (the entry, and the 11-04 weekly tonight). This item supplies that call's exposure list.

**Anchor.** The classification's completeness identity: the per-name totals equal `grep -rnw -c` summed over files, 0 lines are left unclassified, and the file count is 16 at `b29c4f77` (or the HEAD count, with the difference named).
- **Positive control:** `tests/validation/test_imp3_tip_q_entry48_def.py` must classify as **(K), all-cells**, for its gated (P) (item 810). Its retained two-tag residual (`:169`, `:233`) must classify as **(K), two-tag, a control**. That residual is the entry's own symptom, reproduced at 2.8818e-06 (`20261104T104820Z_GEO-36-I57-W1d.log:2010`).
- **Negative control:** one hole-build caller, `tests/validation/test_imp3_rod_birdcage.py`, must classify as **(U)**, with the line that makes its rod a hole cited. If either control comes out otherwise, the rule is wrong: stop and report.

**How.** Code reading only, with no solve and no mesh. For each call, follow the fixture to its build:
- `implant_rod=` / `keep_cells` / `_build_kept_rod` / `_arm_with_403`;
- `BIRDCAGE_IMPLANT_CELL_TAG` in the material map;
- a disc's build, for `test_imp3_disc_birdcage.py` and `examples/implants/09_pec_disc_in_birdcage_rim_sar.py`;
- the `Study` route, for `tests/package/test_wf15_study_implant_rod.py`.

Cite file:line for the deciding statement. A PEC body that is a hole is (U). A body kept as cells at σ > 0 is (K). Where the reading needs the built mesh's tags to decide, write (?) and name the build. **Do not build it.**

**Window.** One smoke window records the count:

`scripts/testing/run_and_log.sh OPS-80-grep "timeout -k 30 30 grep -rnw -e _exact_shares -e _loss_power_w tests src examples scripts"`

It is bare, foreground and repo-relative. Its log is the count's source. No other window is needed.

**Tier.** smoke, no container, no solve. **Class:** ops / docs, ≈ 3.1 measured / 20 predicted min. **Depends on:** nothing. **Family:** `OPS-80`, attempt 1 of 4.

**Traps already paid for.**
- `grep -w` counts `_` as a word character. A wrapper named, say, `_exact_shares_all` is therefore not in the count, while a call through a module attribute (`_pd._loss_power_w`) is. Read each line, never the count alone. (This is why `OPS-78` checked totals per name.) A wrapper found that way is listed and classified too, with the difference from `grep` named.
- An import alias (`from … import _exact_shares as _es`) hides later calls from the grep. Follow every alias and list its calls under the import line.
- The tag a kept body carries is the builder's, not the caller's. `build_four_port_sweep`'s material map has no tag 403 (I41 W1b, `test_imp3_tip_q_entry48_def.py:13–17`). A kept body whose map omits 403 takes the default material, so read what that default's σ is before classifying it (U).
- Do not edit the helper in passing. Other modules import it (the entry's own instruction).

**Scope.** This is an exposure list. **Not claimed:** any re-run, any (P) re-read, a band, a fix, or the weekly's widening decision. The known-issues entry **stays open**. Entries leave only with the commit that fixes them.

**Negative result.** If a (K) two-tag call is found on a gated (P) other than the known one, it is the finding: name it in the entry's Scope with file:line and the log of its last green run, and change nothing. If the controls split, report and stop. No `PREDICTED RED:` line.

**Status it can move.** The known-issues entry's Scope moves from "the callers have not been audited" to the named table, which is either "no other (K) two-tag gate" or the named list. **`OPS-80` ⬜ → ✅** on the completeness identity and both controls.

**Outcome — DONE 2026-11-05 (19:30 CST slot, 2026-11-04 local; `claude-opus-5-5`).** Window `20261105T013028Z_OPS-80-grep.log` (Status 0, 22 s, at `42dd675d`): **46 text lines in 16 files**, the same as at `b29c4f77` (the 18:00 review's commit changed no code). The log also prints 21 `__pycache__` "binary file matches" lines, which are not counted *(corrected from 24 by the 2026-11-05 03:00 review's audit: `:35`, `:36`, `:45–62`, `:97`)*. Per name, by `grep -rnwI` on the same tree: `_exact_shares` 35 lines, `_loss_power_w` 11. No line matches both, so 35 + 11 = 46. **Result: no other (K) two-tag gate.** The only gated two-tag (P) on a kept build is the known one.

*Completeness identity.* Each line is classified exactly once. `_exact_shares`: 1 definition + 8 imports + 10 prose + 16 calls = 35. `_loss_power_w`: 2 definitions + 1 import + 4 prose + 4 calls = 11. In total, 3 definitions + 9 imports + 14 prose + 20 calls = 46, with 0 unclassified.

*The rule's two facts, cited once.* `build_four_port_sweep` gives lossy material to tags 1 and 3 only, and every other tag gets the default σ = 0 (`tests/validation/test_port_birdcage_four_port.py:472–485`). A body is a hole unless `keep_cells` is set: `src/fem_em_solver/io/mesh.py:4940` (`.get("keep_cells", False)`), with 402 then a facet (`:272`). `ROD` and `SPHEROID` carry no `keep_cells` (`tests/mesh/test_birdcage_implant_rod.py:97`, `tests/mesh/test_birdcage_implant_spheroid.py:76–81`). So a call is (K) only if its caller both keeps cells and adds a lossy material for a tag other than 1 and 3. A kept body whose tag is left out of the map is σ = 0, which makes it (U).

| # | Line | Name | Class | Mesh / deciding statement |
|---|---|---|---|---|
| 1 | `tests/validation/test_birdcage_power_identity.py:313` | `_exact_shares` | definition | — |
| 2 | `tests/validation/test_port_drive_superposition.py:205` | `_loss_power_w` | definition | the two-tag sum, `:213–220` |
| 3 | `examples/ports/13_birdcage_asymmetric_drive.py:187` | `_loss_power_w` | definition | a separate local copy (4 args), same two tags |
| 4 | `tests/package/test_wf15_study_implant_rod.py:43` | `_exact_shares` | import | — |
| 5 | `tests/validation/test_birdcage_b1_plus_closed_form.py:195` | `_exact_shares` | import | — |
| 6 | `tests/validation/test_imp3_disc_birdcage.py:110` | `_exact_shares` | import | — |
| 7 | `tests/validation/test_birdcage_power_identity.py:175` | `_loss_power_w` | import | — |
| 8 | `tests/validation/test_geo27_one_gram_rung_degree2.py:55` | `_exact_shares` | import | — |
| 9 | `tests/validation/test_birdcage_power_identity_degree2.py:89` | `_exact_shares` | import | — |
| 10 | `tests/validation/test_ans4_h33_sheet_width_excess.py:224` | `_exact_shares` | import (lazy) | — |
| 11 | `tests/validation/test_imp3_tip_ratio_1g.py:61` | `_exact_shares` | import | — |
| 12 | `tests/validation/test_imp3_rod_birdcage.py:84` | `_exact_shares` | import | — |
| 13 | `tests/validation/test_birdcage_b1_plus_closed_form.py:247` | `_exact_shares` | prose (comment) | — |
| 14 | `tests/validation/test_port_lumped_rlc_termination.py:2182` | `_exact_shares` | prose (comment) | — |
| 15 | `tests/validation/test_birdcage_power_identity_degree2.py:24` | `_exact_shares` | prose (docstring) | — |
| 16 | `tests/validation/test_wf7_f_human_ring_matrix.py:805` | `_loss_power_w` | prose (comment) | — |
| 17 | `tests/validation/test_port_drive_superposition.py:228` | `_loss_power_w` | prose (docstring) | — |
| 18 | `tests/validation/test_port_drive_superposition.py:653` | `_exact_shares` | prose (assertion message) | — |
| 19 | `tests/validation/test_port_drive_superposition.py:1016` | `_exact_shares` | prose (docstring) | — |
| 20 | `tests/validation/test_port_drive_superposition.py:1187` | `_loss_power_w` | prose (docstring) | — |
| 21 | `tests/validation/test_imp3_tip_q_entry48_def.py:27` | `_exact_shares` | prose (docstring) | — |
| 22 | `tests/validation/test_imp3_tip_q_entry48_def.py:199` | `_exact_shares` | prose (docstring of `_two_tag_resid`) | — |
| 23 | `tests/validation/test_imp3_tip_q_entry48_def.py:234` | `_exact_shares` | prose (f-string label) | the value beside it is `_two_tag_resid`, see row 47 |
| 24 | `tests/validation/test_imp3_tip_q_entry48_def.py:375` | `_loss_power_w` | prose (comment) | — |
| 25 | `tests/validation/test_imp3_tip_q_entry48_def.py:468` | `_exact_shares` | prose (f-string label) | the diag arm, see row 47 |
| 26 | `src/fem_em_solver/ports/shares.py:4` | `_exact_shares` | prose (docstring) | — |
| 27 | `tests/validation/test_birdcage_power_identity.py:330` | `_loss_power_w` | call | inside `_exact_shares`, so each `_exact_shares` call below decides it |
| 28 | `tests/package/test_wf15_study_implant_rod.py:256` | `_exact_shares` | call — **(U)** | `_build_rod(…)` hole (`:251`); the Study refuses `keep_cells` (`:191`, `:205`) |
| 29 | `tests/validation/test_birdcage_b1_plus_closed_form.py:583` | `_exact_shares` | call — **(U)** | `build_four_port_sweep(phantom_material=VACUUM, …)` (`:513`), no implant |
| 30 | `tests/validation/test_imp3_disc_birdcage.py:304` | `_exact_shares` | call — **(U)** | `_arm_disc` on `_build(True)` / `_build_spheroid()` (`:485`, `:490`), a disc hole |
| 31 | `tests/validation/test_imp3_disc_birdcage.py:647` | `_exact_shares` | call — **(U)** | `_control_c1_disc` on the same disc-hole case (`:632`) |
| 32 | `tests/validation/test_birdcage_power_identity.py:370` | `_exact_shares` | call — **(U)** | `build_four_port_sweep()` (`:363`/`:365`) |
| 33 | `tests/validation/test_birdcage_power_identity.py:381` | `_exact_shares` | call — **(U)** | the ×2 control on the same mesh (`:378–380`) |
| 34 | `tests/validation/test_geo27_one_gram_rung_degree2.py:210` | `_exact_shares` | call — **(U)** | `_build(True, phantom_resolution=…)` (`:150`), no implant |
| 35 | `tests/validation/test_port_lumped_rlc_termination.py:2323` | `_exact_shares` | call (`port16.`) — **(U)** | `build_four_port_sweep(frequency_hz=…)` (`:2308`) |
| 36 | `tests/validation/test_birdcage_power_identity_degree2.py:168` | `_exact_shares` | call — **(U)** | `build_four_port_sweep(…, build_only=True)` (`:235`) |
| 37 | `tests/validation/test_ans4_h33_sheet_width_excess.py:247` | `_exact_shares` | call — **(U)** | `_four_port_rung` (`:231`), map tags 1/3 only (`test_port_birdcage_leg_offset_sweep.py:380–388`); also reached from `examples/ports/23_lumped_sheet_excess_vs_width.py:145` |
| 38 | `tests/validation/test_imp3_tip_ratio_1g.py:182` | `_exact_shares` | call — **(U)** | `_arm_d` on `_built` (`:121–128`), `_build` or the `_build_rod` hole; also `test_imp3_tip_ladder_1g.py:200`, `:392`, `:620` on `_build_rod` holes (`:196`, `:353`, `:565`) |
| 39 | `tests/validation/test_port_drive_superposition.py:465` | `_exact_shares` | call (`step1.`) — **(U)** | `superposition_case`, `build_four_port_sweep()` (`:314`) |
| 40 | `tests/validation/test_port_drive_superposition.py:482` | `_loss_power_w` | call — **(U)** | same case |
| 41 | `tests/validation/test_port_drive_superposition.py:1374` | `_loss_power_w` | call — **(U)** | `_build_ring_quadrature_case`: the ring context (`test_port_birdcage_ring_column.py:404–410`, tags 1/3) or the WF-7 F-human context (`test_wf7_f_human_ring_matrix.py:756–763`, tags 1/3); callers `test_wf7_f_human_b1_quadrature.py:275`, `test_wf7_f_human_sar10g_hotspot.py:215`, `examples/ports/15_birdcage_sixteen_leg_quadrature_b1.py:170` |
| 42 | `tests/validation/test_imp3_rod_birdcage.py:658` | `_exact_shares` | call — **(U)** | `_control_c1` on the rod arm's case, the `_build_rod` hole (`:515`); also `test_imp3_tip_ratio_1g.py:382` on its hole |
| 43 | `examples/implants/09_pec_disc_in_birdcage_rim_sar.py:287` | `_exact_shares` | call (`g9.`) — **(U)** | `g9._case_disc(None)` (`:213`), the disc hole |
| 44 | `examples/implants/09_pec_disc_in_birdcage_rim_sar.py:312` | `_exact_shares` | call (`g9.`) — **(U)** | the same case with the natural wall |
| 45 | `examples/ports/13_birdcage_asymmetric_drive.py:515` | `_loss_power_w` (local) | call — **(U)** | `build_four_port_sweep()` (`:266`) |
| 46 | `tests/validation/test_imp3_rod_birdcage.py:329` | `_exact_shares` | call, inside `_arm` (`:295`) — **split by caller, rows 46a–46f** | the build and map are the caller's |

The wrappers behind row 46 are not counted by `grep -w`. This is the difference from `grep` the trap names: `_arm` passes the call on under another name.

| Row | Caller of `_rb._arm` | Class | Deciding statement / gate |
|---|---|---|---|
| 46a | `test_imp3_rod_birdcage.py:511`, `:520` (its own W1/W2) | **(U)** — **negative control ✔** | off = `_build(True)`. Rod = `_build_rod(**kw)` (`:515`), whose `ROD` has no `keep_cells` (`test_birdcage_implant_rod.py:97`, `:146`), so the rod is a hole with facet 402 pinned (`:7`, `:22`). (P) is asserted two-tag (`:626–627`), which is correct here |
| 46b | `test_imp3_rod_heating_phantom.py:264`, `:269` (I39) | **(U)** | the `_build(True)` / `_build_rod` holes. `shares` is never read |
| 46c | `test_imp3_tip_q_entry48_def.py:119` via `_arm_rod` (rod and rod2 arms, `:259`, `:263`) and `_gel_seen` (gel twin, `:261`) | **(K), all-cells** — **positive control ✔** | `_build_kept_rod()` (`:252`) with 403 mapped to σ_wire (`:111–113`) or to the gel. The gated (P) is `_p_resid` → `_all_cells_resid` whenever the arm carries per-tag losses (`:168–171`, `:204–215`; gate `:345–347`), green ≤ 8.3e-15 (`20261104T104820Z_GEO-36-I57-W1d.log:1968–1994`) |
| 47 | same module: `_two_tag_resid` (`:198–201`), printed at `:233–235`, asserted by `test_control_two_tag_p1_reproduces_omitted_term` (`:351–357`); diag arm `:449` → `test_diag_p_rod_arm` (`:515–518`, behind `FEM_EM_IMP3_I57_DIAG`) | **(K), two-tag, a control** — **positive control ✔** | the entry's own symptom. The control reproduces 2.882e-06 to 3 digits (`20261104T104820Z_GEO-36-I57-W1d.log:2010`). The diag arm's two-tag (P) is item 801's pre-registered red, the known case, not a new finding |
| 46d | `test_imp3_rod_heating_matched.py:124` via `_arm_with_403(map_gel=True)` (I41 arm M, `:171`) | **(K), two-tag, unread** | `_build_kept_rod()` (`:139`) with 403 mapped to the gel (`:113–116`). `shares` is computed but the module never reads `p_vol`, and asserts and prints no (P) (no `shares` / `p_vol` / `DISCRETE_IDENTITY` in the module). No gate is exposed |
| 46e | `test_imp3_rod_heating_matched.py:176` (`map_gel=False`, the control) and `:173` (arm R) | **(U)** | 403 is unmapped and so σ = 0 (`test_port_birdcage_four_port.py:472`). Arm R is the `_build_rod` hole (`:136`) |
| 46f | `examples/implants/15_rod_heating_matched_phantom.py:164` (`i41._arm_with_403(map_gel=True)`), `:167` (arm R) | **(K), two-tag, unread** / **(U)** | as 46d / 46e. The example reads no `shares` |

Other wrappers followed and found (U): `_exact_shares_w` (`test_port_drive_superposition.py:241`, calls `:462`, `:509`, `:1373`), `_arm_d`, `_control_c1`, `_arm_disc` / `_control_c1_disc`, `_width_case`, `_solve_and_read`, `_build_ring_quadrature_case`. `_exact_shares_enabled` (`test_birdcage_b1_plus_closed_form.py:269`) is an environment knob, not a wrapper. No import aliases the two names: every alias is a module alias (`step1.`, `port16.`, `g9.`, `_rb.`), and each module-alias call is one of the grep rows above.

**Verdict.** The exposed (K) callers are entry48 (gated all-cells, with an asserted two-tag control), I41 arm M and `EX-155` arm M (two-tag shares computed but never read). **No other (K) two-tag gate.** No re-run, no (P) re-read, no band, no edit to the helper or any test. Whether to widen the helper stays the weekly's call. The known-issues entry stays open.
