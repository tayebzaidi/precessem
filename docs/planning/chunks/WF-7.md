# WF-7 — SAR10g hotspot identification

 SAR10g hotspot identification | ✅ (**Step 0d (`xxl`, chain H2) RUN 2026-09-26, recorded by the 03:00 review:** the full degree-2 F-human 32×32 at 64 MHz in one window — one factorisation plus 31 back-substitutions, 3 256 418 unknowns, 10.65 min, 103.06 GiB peak, reciprocity 9.715e-14 asserted; the printed terminal residual is 1.08–1.18e-2 on **every** column (degree 1: 4.1–4.6e-3), unattributed and gating nothing (`20260926T070009Z_WF-7-step0d.log`, `xxl-ledger.md`). The closure below is degree 1 and unaffected. **Audit PASS 2026-09-26** (`20260925T093315Z_WF-7-step3.log:11615–11617`). *Then:* **step 3 DONE 2026-09-25 (09:30 slot, §9 item 39, chain H5): the 10 g SAR hotspot located at human scale, C16- and mirror-gated.** New `tests/validation/test_wf7_f_human_sar10g_hotspot.py` (`FEM_EM_WF7_STEP3=1`) + additive `phantom_resolution=None` keyword on `_build_f_human` (`None` ⇒ no new generator key); no `src/`. `20260925T093315Z_WF-7-step3.log`, Status 0, 352 s, heavy, `-n 8`, 4 passed (`:11666`). Refined phantom 0.0075 m: **537 546 cells, 22 594 phantom cells** (printed beside 504 642, not asserted — no record), build 133.33 s (`:10468`); measured phantom (allreduced tag-3 vertex bbox) centre (0, 7.2e-9, 0) m, R = 64.2857 mm, H/2 = 85.7143 mm; 10 g ball a = 13.3650 mm = 3.56 phantom h across (`:11559`). Nine-ball lattice at az_ref 11.25°, ccw drive: **c\* = r 45 mm, z = z_c**, SAR₁₀g 2.0858e-04 W/kg (ball masses 9.984–10.024 g; the ±50 mm rows peak at 1.44e-04, `:11589–11598`, `:11614`). (i) coverage identity, ball radius 108.1429 mm (< 140 mm conductor clearance, asserted): power rel **9.592e-14**, mass rel **1.297e-13** ≤ `EXACT_IDENTITY_RTOL` 1e-10 (`:11617`); (ii) **C16 worst 0.4904 %** (R_7; range 0.0098–0.4904 %) ≤ `C4_COVARIANCE_BAND` 5 % (`:11599–11615`); (iii) **mirror 0.0502 %** ≤ 5 % (`:11616`) — *disclosed: c\* lies on the az_ref ray, which is the mirror plane, so `_mirror_xy(c*)` = c\* (identical ball mass 10.001116 g, `:11584`); the identity read here is cw SAR at c\* vs ccw SAR at c\*, the pre-registered one, not an off-plane reflection*; (iv) all 30 balls contained 1 mm clear (asserted). Negative control (*predicted* ≫ 5 %, printed, never asserted): single-drive P17 SAR₁₀g at c\* and C4 images R_4/R_8/R_12 68.11 / 96.40 / 70.03 %, spread 96.40 % = 196.6× the ccw worst — as predicted (`:11619`). Printed only: whole-phantom mean SAR (ccw) 1.7356e-04 W/kg (`:11618`); step 2's `|B₁⁺|` on the refined mesh C16 worst 0.0513 %, mirror 0.0456 %, **215 / 215** valid sample points (the 38-point question), identity 3.475e-16 (`:11620`); PRICE build 133.33 s, sweep+case 93.86 s, SAR phase 120.40 s (30 balls, 2.5–4.7 s each — the 25 s STOP rule did not fire), 348.02 s wall, 19.052 GiB summed `ru_maxrss` (`:11621`). Rule (a): the step-2 module re-run with the keyword unset, 5 passed, 507 266 cells, C16 worst 0.1495 % reproduced, 210 s (`20260925T093918Z_WF-7-step3-ruleA-step2.log:10480`, `:11566`, `:11609`). *One phantom, one fixture, 64 MHz, degree 1 — a SAR location and two symmetry identities; no C95.3, compliance or absolute SAR claim; carries the 5.50 % order caveat (step 0b).*) (**step 2, 2026-09-23 (04:30 slot, §9 item 36, chain H4): `|B₁⁺|` C16 gated at human scale.** New `tests/validation/test_wf7_f_human_b1_quadrature.py` (`FEM_EM_WF7_STEP2=1`) + additive `built=None, cell_record=None` keywords and four additive return keys (`result`, `points`, `mask`, `az_ref`, needed by (vi) and the P17 control) on `_build_ring_quadrature_case()`; no `src/`. `20260923T095104Z_WF-7-step2.log`, Status 0, 227 s, `-n 8`, 5 passed. On the F-human fixture at 64 MHz, degree 1, 32-port ccw quadrature with `POST-6` step 3's weights: (i) 507 266 cells (step 1's count reproduced), +0.52 % to the 504 642 record inside `CELL_COUNT_BAND` (`:10480`); (ii) 16 ports per ring on the 22.5° grid (`_ring_quadrature_slots` asserts; reference azimuth 11.250000°, `:11550`); (iii) worst C16 image of `|B₁⁺|_ccw` **0.1495 %** (R_8; range 0.0978–0.1495 %) ≤ `C4_COVARIANCE_BAND` 5 % (`:11551–11566`); (iv) mirror `|B₁⁻|_cw(Mx)` vs `|B₁⁺|_ccw(x)` **0.1634 %** ≤ 5 % (`:11567`); (v) `PORT-16` exact discrete identity on the superposed field 2.316e-16 ≤ 1e-6 (`:11568`); (vi) unit weight e_P17 deviation 0.000e+00 ≤ 1e-12 (`:11571`). Predicted controls, printed only: cw mis-paired 99.2969 % vs predicted 99.3438 % (`:11569`); single-drive P17 worst C16 spread 33.3425 %, 223× the ccw worst (`:11572`). Printed: **38 / 38 sample points — below the 50-point `MIN_SAMPLE_POINTS` floor `POST-6` step 3 asserts on its own fixture, not gated here (review to rule)**; `S[P17,P17]` 0.407423+0.344417j = step 1's digits (`:11573`); PRICE build 134.02 s, sweep 66.20 s, 223.82 s wall, 17.713 GiB summed `ru_maxrss` (`:11574`). Rule (a): `POST-6` step 3 re-run 17 / 17 green, 0.8102 % / 0.6769 % reproduced, 172 s (`20260923T095500Z_POST-6-step3-rerun-WF-7-step2.log:13670–13674`). *C16 invariance of `|B₁⁺|` at human scale only — no homogeneity, absolute `B₁⁺`, SAR or Larmor-accuracy claim.* *Self-consistency only; degree 1 at human scale carries the weekly's 5.50 % order caveat (`WF-7` step 0b, 2026-09-19).* H5 (step 3) becomes queueable. **Ruling on the 38 sample points (2026-09-25 03:00 review): step 2 stands as recorded.** `MIN_SAMPLE_POINTS` = 50 (`tests/validation/test_birdcage_b1_plus_map.py:135`) belongs to `POST-6`'s fixture. Its sample cylinder is absolute (radius 20 mm, half-height 20 mm, `:125–126`) and does not scale with the phantom, so on F-human's 15 mm phantom cells it holds 38 centroids. Step 2's Done-when (this row) names the band, the two identities and the unit weight, not the floor, and the count is disclosed. Not demoted. Step 3 runs on a refined phantom and reprints the same C16 reading on its denser point set as a record.) (**Step 3 (chain H5; §9 item 39, written 2026-09-25 review under rule (4)) — Done-when (added by that review, rule (j)):** on the F-human fixture at 64 MHz, degree 1, 32-port ccw quadrature (step 2's superposition), with the phantom refined to `phantom_resolution` = 0.0075 m (the 10 g ball, a = 13.365 mm at ρ = 1000, ≈ 3.6 phantom cells across — the regime `MAT-4` step 4 gated 10 g at on F-small; at the default 15 mm it is ≈ 1.8, the regime whose 1 g column missed 5 %), all asserted green in one logged window with elapsed: (i) `MAT-4`'s **coverage identity** — one ball enclosing the whole phantom (centre and radius from the allreduced tag-3 bounding box; radius just over the circumradius and clear of the conductor), ρ = 1000 on tag 3 and 0 elsewhere, reproduces `mean_sar(…, subdomain_ids=PHANTOM_CELL_TAG)`'s dissipated power and ρ × the meshed phantom volume within the imported `EXACT_IDENTITY_RTOL` (1e-10); (ii) the **10 g hotspot** — the argmax of `mass_averaged_sar` over a fixed candidate lattice — and its fifteen **C16 images** within the imported `C4_COVARIANCE_BAND` (5 %); (iii) the **mirror identity** — the cw drive's 10 g SAR at the mirrored hotspot vs the ccw drive's at the hotspot — within the same band; (iv) every ball's containment in the phantom. The location table is printed ⇒ **`WF-7` 🟡 → ✅** at the stated grain: *one phantom, one fixture, degree 1, the 5.50 % order caveat; a SAR location and a symmetry identity, not a C95.3 or compliance figure.*) (step 1, 2026-09-21: the F-human degree-1 32×32's identities gated — self-consistency only, see the step-1 record at the end of this cell) Step 0 (cost probe, 2026-09-13): one degree-1 single-drive lumped-sheet solve on the F-human rung (longitudinal ring sheets, 507 266 cells, +0.52 % vs `GEO-25`'s transverse record) at 64 MHz, `-n 8` — 607 039 unknowns, solve 37 s, mesh build 123 s, summed `ru_maxrss` 10.93 GiB (`20260913T190102Z_WF-7-step0.log:10419–10427`, 178 s); below the predicted 11–33 GiB / 3–8 min bracket. A measurement for the XXL commissioning, no physics claim. **Step 0b queued 2026-09-13 (operator, interactive) for the Saturday 2026-09-19 02:00 `xxl` window (`docs/testing/xl-pending.md` entry 1; `docs/testing/xxl-queue.env`):** the same probe with the new `FEM_EM_WF7_DEGREE=2` knob (unset is byte-identical), degree-1 control first in the same window at `-n 16`, against `fem-em-solver-xxl`. Readout, record only: `S_driven` at degree 1 and 2 on the same F-human mesh and `|S₂ − S₁|/|S₁|` — the human-scale order sensitivity — plus the degree-2 price (predicted ≈ 3.2 M unknowns, ≈ 235 GiB, ≈ 17 min from 2d's point; brackets 150–350 GiB / 10–60 min printed). Decision rule for the 09-19 weekly pre-registered in the pending file: same class as F-small's 64 MHz order move ⇒ degree 2 is the production order at human scale and this is its price; an order of magnitude smaller ⇒ degree 1 suffices for human-scale S; does not finish ⇒ the price is the finding and `TH-16` moves up. Nothing asserted beyond the imported cell band. **Step 0c (opened 2026-09-16 review, §9 item 1 + `xl-pending.md` entry 6):** the probe gains `FEM_EM_WF7_PORTS` (unset = byte-identical; `k` or `all` = the first `k` / every ring port driven in turn, each column's fields dropped), proved at heavy tier by the flag-off control (`S_driven` = step 0's `0.407423+0.344417j`) and a two-drive reciprocity assert; then the XL cost probe — the full 32-port degree-1 set on the F-human mesh at `-n 16`, readout the price plus the 32×32's reciprocity / passivity / C16 spreads printed beside `PORT-13`'s bands, nothing asserted. **The knob landed 2026-09-16 (heavy, `-n 8`, two windows): (a) flag-off control — the knob unset, the probe reproduces step 0's printed `S_driven` digits `0.407423+0.344417j` exactly, cells 507 266 (+0.52 %, inside the imported band), solve 31.71 s, summed `ru_maxrss` 10.93 GiB, 158 s window (`20260916T093301Z_WF-7-step0c.log:10425–10430`); (b) `FEM_EM_WF7_PORTS=2` — two drives P17, P18, the 2×2's reciprocity ratio **9.767e-16** against the imported `RECIPROCITY_BAND` 1e-3 and drive 2's `S_driven` 0.407959+0.343115j differing from drive 1's by 1.407e-03 > 1e-6, 155 s window (`20260916T093548Z_WF-7-step0c.log:10431–10438`). Priced for the XL window: the first drive factorises (27.71 s), the second back-substitutes the held MUMPS factor (`solve_kind` `held`, **0.59 s**, `PORT-19` step 3), so 32 drives extrapolate to ≈ 2 min mesh + 28 s + 31 × 0.6 s ≈ 3.5 min at `-n 8`; summed `ru_maxrss` moved 10.892 → 11.082 GiB over the second drive (+0.19 GiB per retained column). Negative control (*predicted*, printed, never asserted): the 2×2's `σ_max` 0.709401 ≤ `COLUMN_PASSIVITY_CEILING` 1 — as predicted; class spreads printed (|S_jj| spread 8.0579e-04, full-column Σ|S_ij|² 0.895804 / 0.895395). No F-human `S` claim, no band moved.** *The 32-port XL cost probe was queued by the 2026-09-18 03:00 review for Monday 2026-09-21 02:00 (`docs/testing/xl-queue.d/20260921-WF-7-step0c.env`, `xl-pending.md` entry 6).* **Step 0b ran 2026-09-19 02:00 (`xxl`, cron, no session) — a record, 🧪, banked by the 03:00 review from a `log-pathologist` reading** (`20260919T070008Z_WF-7-step0b.log`; both legs rc 0, Status 0, 673 s, `memory.peak` 106.1 GiB): degree 1 at `-n 16` reproduces step 0's `S_driven` digits (printed, not asserted at this width; 607 039 unknowns, 24.61 s); **degree 2 is 3 256 418 unknowns, 347.56 s, 107.9 GiB summed** — a third of the predicted price on both axes, so a degree-2 F-human drive is `xl`-sized; `S_driven(P17)` `0.436656+0.346994j`, i.e. `|S₂ − S₁|/|S₁|` = **5.50 %** (the review's arithmetic from two printed values; nothing about the degree-2 `S` is asserted). Unattributed and for the weekly: the printed power-accounting residual is 1.104e-02 at degree 2 against 4.187e-03 at degree 1. The decision rule's branch (c) is excluded; (a)-vs-(b) is the 2026-09-19 weekly's. **Step 0d** (daily licence, priced-family variant: degree 2 at the full 32-port set) is queued for the 2026-09-26 `xxl` window (`xl-pending.md` entry 9). Probe label defects from this log: `OPS-52`. **Step 0c RUN 2026-09-21 (`xl`, chain H1; `20260921T070008Z_WF-7-step0c.log`, Status 0, 228 s, `log-pathologist` CONFIRMED countable):** the degree-1 32-drive set is one factorisation (24.73 s) + 31 held back-substitutions (0.40–0.57 s), 607 039 unknowns, cumulative 224.19 s, summed `ru_maxrss` 13.842 GiB, `memory.peak` 11.96 GiB at `-n 16`; full 32×32 reciprocity **1.084e-14 ≤ the imported 1e-3, asserted**; σ_max 0.999915, `|S_jj|` spread 1.7891e-03, column sums 0.895287–0.896108, terminal residual 4.061e-03–4.624e-03 — printed only. ⇒ the F-human degree-1 32×32 is a **heavy-tier** job. **Step 1 (chain H3; §9 item 21, written 2026-09-21 review under rule (4)) — Done-when (added by that review, rule (j)):** on the F-human fixture at 64 MHz, degree 1, all 32 drives: the cell band, reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ `MATRIX_PASSIVITY_CEILING`, exactly 18 C16 × mirror classes with worst spread ≤ `OPPOSITE_SPREAD_BAND`, `PORT-16`'s exact discrete power identity ≤ `DISCRETE_IDENTITY_RTOL` on the driven columns, and the mesh-level mis-paired-orientation control, all imported and asserted green in one logged window with elapsed ⇒ `WF-7` 🧪 → 🟡. *Self-consistency only; degree 1 at human scale carries the weekly's 5.50 % order caveat (`WF-7` step 0b, 2026-09-19).* **Step 1 DONE 2026-09-21 (09:30 slot, heavy, `-n 8`, one window 372 s, Status 0; new `tests/validation/test_wf7_f_human_ring_matrix.py`, no `src/`; `20260921T093429Z_WF-7-step1.log`):** (i) 507 266 cells, +5.200e-03 vs 504 642 inside `CELL_COUNT_BAND` 0.01 (`:10480`); 1 factorisation (33.65 s) + 31 held (0.51–0.70 s), summed `ru_maxrss` 11.012 GiB (`:10484–10516`); (ii) reciprocity **9.200731e-15** ≤ 1e-3 (`:10518`); (iii) σ_max **0.999914648** ≤ `MATRIX_PASSIVITY_CEILING` 1.0, margin +8.535e-05 (`:10519`); (iv) **18** C16 × mirror classes from measured azimuths / sheet-centre z, worst azimuth residual 5.684e-14°, worst spread **3.6997 %** (other ring, 2 steps) ≤ 5 % — first measurement on this fixture (`:10520–10539`); (v) exact discrete power identity max **2.835e-15** ≤ 1e-6 on **drives 1, 9, 17, 25 only** — the first identity evaluation cost 17.31 s > the item's 15 s stop rule, so the item's own fallback ran (`:10481–10483`, `:10540`); P_src,exact > 0 on each. Controls: *asserted* default-orientation build (504 642 cells, 132.90 s, no solve) max sheet extent along φ̂ 2.816e-17 m < 1e-12 while the longitudinal chord / min extent is 1.716154e-02 m (`:20224`); *asserted* P17 column × 1.01 lifts the ratio to 2.064468e-03 = 2.244e+11× the measured (bar 1e6×, `:10542`); *predicted, printed* 2.064468e-03 vs the 2e-3 `MATRIX_CONTROL_MARGIN` bar — at/above, inside the review's 1.5–2.1e-3 prediction (`:10543`). Printed only: terminal residual 4.061e-03–4.624e-03 beside pooled C/terminal−1 4.591e-03–4.689e-03; column sums 0.895287–0.896108; S_driven(P17) 0.407423+0.344417j = the step-0 record digits at `-n 8` (`:10541`). ⇒ `WF-7` 🧪 → 🟡. **Self-consistency identities on one fixture at degree 1 — no field, B₁⁺, SAR or absolute S; degree 1 at human scale carries the weekly's 5.50 % order caveat (`WF-7` step 0b, 2026-09-19).** **Step 2 (chain H4; §9 item 36, written 2026-09-23 review under rule (4)) — Done-when (added by that review, rule (j)):** on the same F-human fixture at 64 MHz, degree 1, the 32-port ccw quadrature drive superposed through `superpose_drives` with `POST-6` step 3's weights: the worst of the fifteen C16 images of the CG1 `|B₁⁺|_ccw` map, the mirror identity `|B₁⁻|_cw(Mx)` vs `|B₁⁺|_ccw(x)` — both ≤ the imported `C4_COVARIANCE_BAND` — `PORT-16`'s exact discrete identity on the superposed field ≤ `DISCRETE_IDENTITY_RTOL`, and the unit-weight identity (`UNIT_WEIGHT_CONTROL_RTOL`), asserted green in one logged window with elapsed, plus the rule-(a) re-run of `POST-6` step 3 green ⇒ the row records step 2; `WF-7` stays 🟡 (step 3, the 10 g hotspot, is the chunk's title). *C16 invariance of `|B₁⁺|` at human scale only — no homogeneity, absolute `B₁⁺`, SAR or Larmor-accuracy claim; carries the 5.50 % order caveat.* 

### §9 item 213 — 2026-10-04 (chain H15a)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.*

**`WF-7` step 0e, prerequisite half (§10 H15a): a save knob `FEM_EM_WF7_SAVE_NPZ` on `tests/validation/test_wf7_f_human_ring_matrix.py`, with the degree-1 32×32 reproduced through the saved file at the imported reproduction band.** It unlocks `docs/testing/xl-pending.md` entry 29 (H15b, degree 2, κ off, saved).
- **Executor:** implementer. **Tier:** heavy, `-n 8`, complex (W1); W0 standard, `-n 1`. **Class:** chain step (solve), measured median ≈ 18 min.
- **Cost:** ≈ 18 measured / ≈ 40 predicted min (15 fixed + ≈ 15 code + W0 ≈ 1 min + W1 ≈ 7.5 min + ≈ 2 spare).
- **Depends on:** nothing. The module is as `POST-7` step 1 left it (knobs at `:168–187`); check `git log -1 -- tests/validation/test_wf7_f_human_ring_matrix.py` before editing.
- **Family:** `WF-7` step 0e (H15), attempt 1 of 4.

**The contract with entry 29 (binding).** The knob is named exactly `FEM_EM_WF7_SAVE_NPZ`, beside `FEM_EM_WF7_DEGREE` and `FEM_EM_WF7_IDENTITY_COLUMNS`. The implementer who lands this green edits **only entry 29's Status line**, in the same commit: PENDING PREREQUISITE → READY, citing W1's log and the knob's value convention (below). It does **not** write, issue, queue or edit entry 29's command, and it does not touch `docs/testing/xl-queue.d/`. The command is the review's to issue ("entry 21's shape with `FEM_EM_WF7_DEGREE=2`, all 32 columns, the save knob set"). If any other name lands, the review re-issues entry 29; nobody edits it.

**Code (tests only, no `src/`; additive; rule (a), with W1 as the re-gate).** All changes go in `test_wf7_f_human_ring_matrix.py`.
- **`_save_knob()`** is a new function. Do not change `_post7_knobs()`'s return signature. It reads the env at call time:
  - unset or empty ⇒ `None`: no file, and behaviour otherwise unchanged;
  - otherwise the value is an **absolute path under `/workspace/logs/` ending `.npz`** (e.g. `/workspace/logs/wf7/wf7-s32-degree1-n8.npz`). Anything else raises `ValueError` before the mesh is built;
  - set together with `FEM_EM_WF7_IDENTITY_COLUMNS` ⇒ `ValueError`. The file holds only the full 32×32.
  - Add the value to the `start:` line's knob print.
- **Save point:** inside `if matrix:`, straight after `s` is formed (`:412–415`). That is before any further print or assert, before the default-orientation control build (`:476–496`, ≈ 125 s), and before the asserts (`:500–550`). So a degree-2 red at (ii)–(iv), or a kill during the control build, still leaves the matrix on disk.
- **Write:** rank 0 runs `os.makedirs(dirname, exist_ok=True)` and `np.savez(...)`, then `comm.bcast`s an ok flag; all ranks hit `comm.Barrier()`. Every rank raises if the write failed, so no rank is left waiting in a Barrier.
- **Keys:** `ans:8`'s schema (`examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py:241–245`), so `_slot_port_map` / `ring_mode_rows_v2` read it unchanged:
  - `s32` (complex128, rows and columns in `order`), `order` (`driven_ids`, P17…P48, a numpy unicode array);
  - `slot_k`, `slot_sign` from `_ring_quadrature_slots` (`test_port_drive_superposition.py:1208`) on `[{"ordinal": i, "tag": SHEET_IFACE + i, "azimuth_deg": azimuth[i], "z": centre_z[i]} for i in ring_ports]`;
  - `kappa` (P17's pooled C/terminal − 1, a diagnostic), `cells`, `freq_hz`, `degree`, `reciprocity`, `sigma_max`;
  - plus `kappa_applied=False` (the module's widths are already area/chord, `:255`), `pooled_kappa` per column (NaN where the identity was not evaluated), and `ranks`.
- **Reload, "through it":** after the Barrier, every rank runs `np.load(path, allow_pickle=False)`.
  - The class means are computed from the **reloaded** array by a new `_class_means(s, az, zc)`, lifted verbatim from `:423–439`.
  - The inline (iv) block is **not** rewritten, so the gate code path stays as it was.
- **W0's test:** a new `test_wf7_save_knob_contract`. It uses `monkeypatch` env only, builds no mesh, and is not `complex_only`.
- **Leave alone:** `_sheet_centre_z`. `test_wf7_f_human_b1_quadrature.py:112` imports it, so every `WF-7` step 2/3 and `WF-12` module depends on it.

**Anchors (asserted).**
- Every existing anchor stays unchanged and is re-gated by W1 under rule (a): (i)–(v), the P17 reproduction at `REUSE_REPRODUCTION_RTOL`, and both asserted controls.
- **(L1) Round trip:** `s32.dtype == complex128`, `np.array_equal(loaded, s)`, `order` equal, and 16 slots per ring in the stored map.
- **(L2) Degree-1 reproduction through the file:** from the reloaded array, σ_max against `SIGMA_MAX_D1_RECORD = 0.999914648` (`20260928T093759Z_POST-7-step1-W1.log:10641`; equal to step 1's `20260921T093429Z_WF-7-step1.log:10519`), and each of the 18 C16 × mirror class means \|S\| against `CLASS_MEAN_D1_RECORD` (18 literals, `…POST-7-step1-W1.log:10643–10660`, entered with their lines).
  - Each must be within the imported `REUSE_REPRODUCTION_RTOL` 1e-6 (`test_port_birdcage_ring_matrix.py:171`), asserted at degree 1 and printed at degree ≠ 1, as the module already does for P17 (`:394–404`).
  - The records' print resolution is ≤ 5e-10 relative, 2000× under the band. The fixture is deterministic at `-n 8`: P17 reproduced across two windows to 2.4e-11 (`…W1.log:10638`, `…W2.log:10606`).
- **W0:** `_save_knob()` is `None` when unset. It raises on a relative path, on a path outside `/workspace/logs/`, on a non-`.npz` suffix, and when `FEM_EM_WF7_IDENTITY_COLUMNS` is set. No file is written.

**Negative controls (rule (e)).**
- *Asserted (arithmetic, on the reloaded array):* scaling column P17 by `CONTROL_COLUMN_SCALE` (1.01) must move the worst class-mean reproduction by ≥ 10 × `REUSE_REPRODUCTION_RTOL`.
  - **Ceiling:** one column supplies 1 of the 32 entries in an n = 32 class, or 2 of the 64 in an n = 64 class, so the shift is ≈ 0.01/32 ≈ 3.1e-4, i.e. ≈ 310× the band.
  - The 10× bar sits 31× under that ceiling. Predicted range [2e-4, 4e-4].
- *Predicted, printed, never asserted:* the same reproduction on `s.astype(np.complex64)` moves the worst class mean by ≈ 1e-8 to 1e-7. That is **below** the 1e-6 band. The band alone cannot see a single-precision save, which is why (L1) asserts the dtype and bit-equality.
- *Binary:* W0 shows that the unset knob writes no file.

**Printed against predictions written first.** W1 should reproduce the 09-28 W1 digits:
- cells 507 266;
- reciprocity ≈ 1.29e-14, σ_max 0.999914648, worst spread 3.6997 %;
- (v) ≤ 2.6e-15 on all 32 drives (the stop rule does not fire: measured projection 389 s ≤ 540 s), NEG ≈ 479×;
- npz ≈ 17 KiB, save plus reload < 1 s.

**Tier and windows.** Bash `timeout: 660000`, foreground. Each harness call is bare and whole. Run with `-s`.
- **Measured basis, `POST-7` step 1 W1** (this module, knobs unset, degree 1, `-n 8`, all 32 drives with the identity):
  - 437 s wall, pytest 434.8 s, window 399.8 s, summed `ru_maxrss` 11.254 GiB;
  - build 124.16 s, factorising drive 29.67 s, held drives 0.50–0.66 s, identity 2.10 s per drive, control build ≈ 125 s (`20260928T093759Z_POST-7-step1-W1.log:10600–10664, 20346, 20870`).
  - W2 (one column): 279.7 s window, 316 s wall, 11.03 GiB (`20260928T094524Z_POST-7-step1-W2.log:20289, 20813`).
  - `WF-7` step 1 (four-drive identity fallback): 372 s, 11.012 GiB.
- **W0**, standard, `-n 1`, predicted ≤ 60 s:
  `scripts/testing/run_and_log.sh WF-7-H15a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_save_knob_contract -v -s --tb=short'"`
- **W1**, heavy, `-n 8`, predicted [420, 480] s wall and ≈ 11.3 GiB. It uses §5.1 durable capture because it runs within ≈ 2 min of the 590 s ceiling:
  `scripts/testing/run_and_log.sh WF-7-H15a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf7 && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_SAVE_NPZ=/workspace/logs/wf7/wf7-s32-degree1-n8.npz timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short > /workspace/logs/wf7-h15a-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/wf7-h15a-w1-raw.log; cat /workspace/logs/wf7-h15a-w1-raw.log; exit \$rc'"`
  - W1 also collects W0's test, which is cheap, so `monkeypatch` must not leak env into the heavy test.

**Traps.**
- **Save before anything can raise** (see "Save point"). Rank 0 writes, then bcast, then Barrier.
- Load with `allow_pickle=False`; `order` is a numpy unicode array, as `ans:8` stores it.
- `/workspace/logs/` is the gitignored top-level `/logs/` (`.gitignore:82`). Never write under the tracked `docs/testing/logs/`.
- The knob never bypasses the identity stop rule (`:347–360`).
- `-s`, or the prints are captured. Copy the durable-capture idiom verbatim, trailing `; exit $rc` kept. No pipe or `| tail` in the harness call.
- Status 124 is the finding. Run `pgrep -c python3` in the container, then the cache sweep. Never lengthen the timeout.
- No module edit after W1 starts (rule (i)).

**Scope.** Degree 1, `-n 8`, σ = 800 S/m coil, 64 MHz, the existing fixture. **Not covered:**
- degree 2 (H15b, entry 29);
- `C_mode` or any `ring_mode_rows_v2` readout;
- `ANS-8` rows, copper, any band;
- issuing entry 29's command.

**Status it can move.**
- `xl-pending.md` entry 29: PENDING PREREQUISITE → **READY** (Status line only, same commit).
- The §7 `WF-7` row records the knob (✅ unchanged).
- §10 H15a lands, so H15b becomes the review's to queue.

**Negative result.**
- **(L1) red:** a save-path defect. Park on `attempt/`, mark BLOCKED (rule (d)); entry 29 stays PENDING.
- **(L2) red with (L1) green:** the degree-1 fixture has moved since 09-28, through an image or mesher change. Stop and report all 18 deltas, BLOCKED. Re-recording falls under the (1*) licence (`record-reconciler`) and never happens in-slot; entry 29 stays PENDING.
- **Any existing anchor red:** the module moved. Stop and park.
- **W1 over 540 s, or Status 124:** the price is the finding. Report it.
- There is no `PREDICTED RED` line.

**Outcome (2026-10-04 04:30 slot, take-next, implementer):** BLOCKED, parked on `attempt/WF-7-step0e-20261004T1007Z` (test module only). W0 `20261004T095522Z_WF-7-H15a-W0.log` Status 0, 13 passed, 100 s. W1 `20261004T095709Z_WF-7-H15a-W1.log` **Status 124 at 591 s**, no Python traceback (PETSc signal 15 at `:11552`, during the orientation-control build) — an overrun: the box ran ≈ 2× slow against 09-28 (build 239.36 s vs 124.16 s, factorising drive 48.24 s vs 29.67 s, first identity 7.69 s vs 2.10 s), so the identity stop rule fired (projection 706 s > 540 s; (v) on drives 1 / 9 / 17 / 25 only). Printed before the kill, asserts not reached: cells 507 266, reciprocity 8.06e-15, σ_max 0.999914648, worst spread 3.6997 %, (v) max 1.97e-15, NEG 480.7×, P17 reproduction 2.4e-11; npz 20 823 B, save + reload 0.027 s; (L1) all True; (L2) σ_max 2.6e-10, worst class mean 4.3e-10 (band 1e-6); column-scale control 3.147e-4 = 314.7× the band (predicted [2e-4, 4e-4]); float32 control 1.27e-8 (predicted 1e-8 to 1e-7). `pgrep -c python3` 0 after the kill. Entry 29 stays PENDING PREREQUISITE. Next: re-run W1 unchanged from the branch when the box is quiet (the 09-28 price was 437 s), or trim the default-orientation control build out of this window; never a longer timeout.

**For the review (outside this item).** W1's npz is a **κ-off** degree-1 32×32. Entry 29's readout puts degree 2 beside `20261002T050705Z_ANS-8-v2-W1-R.log`'s degree-1 rows. Those were built κ-corrected (`ans:8` stage S, `08_…py:221–222`), so read that way the order move and κ are confounded, the same caveat item 175 printed. Running `ring_mode_rows_v2` on this npz gives the κ-off degree-1 comparand at no solve cost. Whether to do that is the review's decision.

### §9 item 227 — 2026-10-04 (chain H15a, re-pose of item 213, attempt 2 of 4)

**Outcome (2026-10-04 12:00 slot, implementer):** BLOCKED, parked on `attempt/WF-7-step0e-r2-20261004T1714Z` (module only, commit `57789717`). W0 `20261004T170126Z_WF-7-H15a-r2-W0.log` Status 0, **14 passed**, 60 s (as predicted). W1 `20261004T170232Z_WF-7-H15a-r2-W1-orientation.log` **Status 124 at 590 s — an implementer defect, not price:** 7 of 8 ranks printed `3 passed, 1 skipped in 131.6 s` (`:9818–9836`; the matrix test skipped before any build, as `only` requires); rank 0 never printed its `orientation control:` line because `test_wf7_default_orientation_control` calls `_rss_summed_gib(comm)` (an `allreduce`) inside `if comm.rank == 0:`, so rank 0 waited in a collective the other ranks never entered. The default build itself completed: one `[birdcage-mesh] mesh wall time 120.57 s` (`:9794`), and pytest's per-rank total 131.6 s ⇒ **B ≈ 121–132 s, a normal-speed box** (09-21: 132.90 s), so the price gate would have passed (≤ 285 s); the extent value and its assert were not printed on rank 0 (the other ranks passed the assert, but the reading is not on the log). `pgrep -c python3` 0 after the kill; cache sweep `20261004T171254Z_WF-7-H15a-r2-cache-sweep.log` found no orphans. W2 not run (rule (i): no module edit after W1 starts; the stop rule allows no third heavy window). Entry 29 stays PENDING PREREQUISITE. **Next (attempt 3 of 4):** on the branch, compute `rss = _rss_summed_gib(comm)` on every rank before the rank-0 print, then run W0 / W1 / W2 exactly as written above; the box was normal-speed at 12:00 CDT.

*Drafted by a general-purpose agent for the 2026-10-04 10:30 daily review; the review reviews it.*

**`WF-7` step 0e, prerequisite half, re-posed (§10 H15a): the same target as item 213, with the default-orientation control build moved into its own window by an additive knob, so neither window can overrun at the ≈ 2× slowdown measured on 10-04.** The target is unchanged: `FEM_EM_WF7_SAVE_NPZ` on `tests/validation/test_wf7_f_human_ring_matrix.py`; the reloaded degree-1 32×32 reproduces σ_max and the 18 C16 × mirror class means within the imported `REUSE_REPRODUCTION_RTOL`; `xl` entry 29 goes READY.
- **Executor:** implementer. **Tier:** heavy, `-n 8`, complex (W1, W2); W0 standard, `-n 1`. **Class:** chain step (solve), measured median ≈ 18 min.
- **Cost:** ≈ 18 measured / ≈ 45 predicted min (15 fixed + ≈ 12 code + W0 ≈ 2 + W1 ≈ 3–5 + W2 ≈ 5–8 + ≈ 3 spare).
- **Depends on:** nothing. Main has not touched the module since the branch: `git log -1 main -- tests/validation/test_wf7_f_human_ring_matrix.py` is still `bfa8232c` (`POST-7` step 1), and the branch's merge base is `0054fe9d`. Re-check before starting.
- **Family:** `WF-7` step 0 (0, 0b, 0c, 0d, 0e). The last sub-step that landed a gate is 0d: a registered record, `xxl` 2026-09-26, with reciprocity asserted. Since then: item 213 (BLOCKED, Status 124) = 1; **this item = attempt 2 of 4.**

**Why item 213 died, measured** (`20261004T095709Z_WF-7-H15a-W1.log`). Every save-path reading printed in band. The window was killed at 591 s by price alone:
- Longitudinal build 239.36 s (`:10617`) against 124.16 s on 09-28.
- Factorising drive 48.24 s, first identity 7.69 s (`:10618, :10621`). So the stop rule fired (projection 706 s > 540 s), and (v) ran on drives 1 / 9 / 17 / 25 only (`:10619–10620`).
- Held drives 1.04–1.29 s, but ≈ 3.9 s of wall per drive (the `t=` stamps `:10622–10652`). Matrix complete and npz saved at t ≈ 447 s (`:10652–10657`).
- The default-orientation control build then started. It reached gmsh 1-D meshing and was killed there (`:11550–11552`), with no assert run. That build is 132.90 s (09-21, `20260921T093429Z_WF-7-step1.log:20224`) / ≈ 125 s (09-28) on a normal box; at 2× it is ≈ 260 s.
- `tests/environment` + collection cost ≈ 35 s on 09-28 (pytest 434.8 s vs window 399.8 s); that is ≈ 70 s at 2×.
- Sum at 2×: ≈ 70 + 447 + 260 ≈ 780 s. No timeout fixes that, and `-n 12` does not help: the build is gmsh-serial, and `ANS-4` item 212 W2 measured `-n 12` *doubling* mesh time on this box. So the control comes out of the matrix window.

**Start from the branch.** `git checkout attempt/WF-7-step0e-20261004T1007Z -- tests/validation/test_wf7_f_human_ring_matrix.py`. This is a path checkout of code, not of a record. `git diff --stat` must show that one file. The save code, (L1), (L2), the column-scale and complex64 controls and `test_wf7_save_knob_contract` are kept **as parked**. Item 213's brief binds them: entry-29 contract, save point, keys, reload, `_class_means`, traps.

**Additive code (tests only, no `src/`).**
- **`_orientation_knob()`** reads `FEM_EM_WF7_ORIENTATION_CONTROL` at call time:
  - unset or empty ⇒ `"inline"`. The landed behaviour: the control is built inside the matrix test after the matrix, as today.
  - `split` ⇒ the matrix test does **not** build the default mesh and does not run the one assert that reads it (`max(default_extent) < DEGENERATE_EXTENT_M`, branch `:749–752`). It prints `orientation control: SPLIT — asserted by test_wf7_default_orientation_control in its own window`. It still asserts the two longitudinal checks (`chord`, `min(long_extent)` > `DEGENERATE_EXTENT_M`, branch `:747–748`), because those read the longitudinal mesh it already holds.
  - `only` ⇒ the matrix test calls `pytest.skip` before any mesh build, and the new standalone test runs.
  - Anything else raises `ValueError`. So does `only` together with `FEM_EM_WF7_SAVE_NPZ` (there is no matrix to save).
- **`_default_orientation_extent(comm, ring_ports)`** is lifted verbatim from branch `:684–698`. It contains the default `_params(…)` with its `"ring_sheet_orientation" not in` assert, `birdcage_port_domain`, `_interface_facet_tags`, and `_sheet_extent_along` against `phi_hat` from `_ring_gap_frame(i, LEG_COUNT)` (pure geometry, the same call as `:366`). It returns `(default_extent, default_cells, build_s)`. **Both** the inline path and the standalone test call it, so W1 executes the very lines the inline path runs.
- **`test_wf7_default_orientation_control`** (`@complex_only`). It skips unless the knob is `only`. It calls the helper, prints the same `orientation control:` line as `:699–703` plus `build_s`, and asserts `max(default_extent) < DEGENERATE_EXTENT_M`.
- **The identity stop rule is not edited.** Under `split` its projection still adds `CONTROL_BUILD_ESTIMATE_S` (135 s), which overstates the window. The fallback therefore fires earlier, never later. Disclosed, and conservative.
- **W0 contract:** a new `test_wf7_orientation_knob_contract` (the name is fixed, because W0 calls it by node id), with monkeypatch only and no mesh. Unset/empty → `"inline"`; `split` / `only` accepted; `"yes"` / `"0"` raise; `only` + save path raises. `test_wf7_save_knob_contract` stays as parked.
- **Leave alone:** `_sheet_centre_z` (imported by `test_wf7_f_human_b1_quadrature.py:112`) and every knob's unset path.

**Anchors (asserted).** Every anchor of item 213 is unchanged:
- (i)–(v);
- the P17 reproduction at `REUSE_REPRODUCTION_RTOL`;
- the column-scale reciprocity control;
- **(L1)** round trip: dtype, `array_equal`, `order`, 16 slots per ring;
- **(L2)** σ_max against `SIGMA_MAX_D1_RECORD` and the 18 `CLASS_MEAN_D1_RECORD` class means, each within the imported `REUSE_REPRODUCTION_RTOL` 1e-6, through the reloaded file.

The mesh-level orientation control keeps its predicate. Only its window moves: W1 asserts it, and W2 asserts its longitudinal half.

**Negative controls (rule (e)).**
- *Asserted, arithmetic, on the reloaded array:* P17 × 1.01 moves the worst class-mean reproduction by ≥ 10 × the band. The ceiling is 0.01/32 ≈ 3.1e-4 ≈ 310×. 10-04 measured **3.147e-4 = 314.7×** (`:10678`). Predicted range [2e-4, 4e-4].
- *Asserted, mesh-level (W1):* the default-orientation sheets have extent along φ̂ < 1e-12 m. 09-21 measured 2.816e-17 m against the longitudinal minimum 1.716154e-02 m, a separation of ≈ 6e14.
- *Predicted, printed, never asserted:* the complex64 copy moves the worst class mean by ≈ 1e-8 to 1e-7, below the band. 10-04 measured 1.271e-8 (`:10679`). This is why (L1) asserts dtype and bit-equality.
- *Binary (executor reads W2's log):* exactly **one** `[birdcage-mesh] mesh wall time` line appears in W2. A second one means `split` did not reach the code.

**Printed against predictions written first.**
- **W1:** default cells 504 642. Extent along φ̂ ≈ 3e-17 m. Build 125–133 s on a normal box; ≈ 2× that if 10-04's load persists.
- **W2** reproduces 10-04's digits:
  - cells 507 266, reciprocity 8e-15 to 1.3e-14, σ_max 0.999914648, worst spread 3.6997 %;
  - (v) max ≤ 2.6e-15 on all 32 drives on a normal box, or on drives 1 / 9 / 17 / 25 if the stop rule fires (both are the module's landed behaviour);
  - NEG ≈ 480×; (L2) worst ≤ 5e-10 relative; npz 20 823 B; save plus reload < 0.1 s.

**Windows.** Bash `timeout: 660000`, foreground. Every harness call bare and whole, `-s` everywhere. **No module edit after W1 starts** (rule (i)).
- **W0**, standard, `-n 1`. Measured 100 s on 10-04 (`20261004T095522Z_WF-7-H15a-W0.log`, 13 passed):
  `scripts/testing/run_and_log.sh WF-7-H15a-r2-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_save_knob_contract tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_orientation_knob_contract -v -s --tb=short'"`
  - Node ids, not `-k`: a `-k` filter would also deselect `tests/environment`.
  - Predicted 14 passed.
- **W1**, heavy, `-n 8`: the control alone **and the box-speed probe**. Predicted 140–150 s on a normal box, ≈ 280 s at 2×:
  `scripts/testing/run_and_log.sh WF-7-H15a-r2-W1-orientation "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_ORIENTATION_CONTROL=only timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short'"`
- **Price gate (after W1, binding).** Read W1's printed default-build seconds, B.
  - W2's projected wall is ≈ 15 + (B / 129) × 450 s. This scales the 10-04 matrix-complete time (≈ 447 s at B-equivalent ≈ 2 × 129) linearly, plus startup without `tests/environment`.
  - **B ≤ 285 s** (slowdown ≤ 2.2×, projection ≤ ≈ 510 s) ⇒ run W2.
  - **B > 285 s** ⇒ **stop**: do not start W2. Commit nothing to the module on `main`. Park on `attempt/` with W0 / W1 logs and the measured B; BLOCKED, "price — box > 2.2× slow". Entry 29 stays PENDING.
- **W2**, heavy, `-n 8`: the matrix, the save, (L1) / (L2), control split out, no `tests/environment` (W0 ran it in this slot). §5.1 durable capture, copied verbatim with the trailing `; exit $rc`. Predicted 290–330 s on a normal box and [440, 520] s at 2×:
  `scripts/testing/run_and_log.sh WF-7-H15a-r2-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf7 && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_ORIENTATION_CONTROL=split FEM_EM_WF7_SAVE_NPZ=/workspace/logs/wf7/wf7-s32-degree1-n8.npz timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short > /workspace/logs/wf7-h15a-r2-w2-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/wf7-h15a-r2-w2-raw.log; cat /workspace/logs/wf7-h15a-r2-w2-raw.log; exit \$rc'"`
  - W2 overwrites 10-04's npz at the same path. That is expected: (L1) and (L2) read the fresh file.
- **Rule (a).** W1 + W2 re-run every test in the edited module under the two knobs, and between them they execute every original line. The unset (`inline`) path is the original control lines plus a call into the lifted helper, and it is **not** run in-slot: an unset degree-1 window is ≈ 440 s normal and ≈ 780 s at 2×. Its first run is entry 29's `xl` window, which asserts the inline control at degree 2. Disclosed; the review rules.

**Stop rule.** One W2 attempt per slot.
- Status 124 on any window is the finding: run `pgrep -c python3` in the container first, then the cache sweep (rubric item 4, only when the count is 0). Park, BLOCKED, with the measured build and drive times.
- Never a longer timeout, never `-n 12`, never a third heavy window in this slot.
- **A third consecutive price death is not a fourth attempt by default.** If this item dies on price, the next review must either queue W2's command as an `xl` priced-family variant of entry 21 (the same module, degree 1, `-n 8`), or wait for a quiet box. It does not reletter.

**Traps** (item 213's, plus):
- the knob must not move the stop rule (`:463–479`) or the `del fields` memory discipline;
- `only` must skip **before** `birdcage_port_domain`, or W1 pays both builds;
- W1 and W2 also collect both contract tests with the knobs set in the environment, so `test_wf7_orientation_knob_contract` must `monkeypatch.delenv` both knobs first, as the parked save contract already does;
- the helper's `assert "ring_sheet_orientation" not in kwargs_default` travels with it;
- W2 without `tests/environment` still sets `FEM_EM_REQUIRE_COMPLEX=1`, so a real-mode container fails loudly rather than skipping;
- rank 0 writes, then bcast, then Barrier (parked code, unchanged);
- never write under the tracked `docs/testing/logs/`;
- no pipe or `| tail` in any harness call.

**Scope.** Degree 1, `-n 8`, σ = 800 S/m coil, 64 MHz, the existing fixture. **Not covered:** degree 2 (H15b, entry 29); `C_mode` / `ring_mode_rows_v2`; `ANS-8` rows; copper; any band; issuing or editing entry 29's command; the unset-path run (see rule (a)).

**Status it can move.**
- `docs/testing/xl-pending.md` entry 29: **PENDING PREREQUISITE → READY**. Status line only, in the same commit. It cites W2's log and the knob's value convention, and adds one clause: `FEM_EM_WF7_ORIENTATION_CONTROL` unset is the landed inline path.
- The §7 `WF-7` row records the two knobs (✅ unchanged).
- §10 H15a lands, so H15b becomes the review's to queue.
- *For the review, not the slot:* entry 29 may want `FEM_EM_WF7_ORIENTATION_CONTROL=split` to save ≈ 2–4 min of `xl` wall. That is a command change, and the review's to make.

**Negative result.**
- **W0 red:** a knob-contract defect. Fix once in-slot before W1; if still red, park.
- **W1 red** (the default build has φ̂ extent): the generator moved. Stop, park, open a known-issues entry; W2 does not run.
- **Price gate trips:** see "Price gate". BLOCKED on price, attempt 2 spent.
- **(L1) red:** a save-path defect. Park, BLOCKED; entry 29 stays PENDING.
- **(L2) red with (L1) green:** the degree-1 fixture moved since 09-28. Report all 18 deltas, BLOCKED. Re-recording is `record-reconciler`'s under (1*), never in-slot.
- **Any (i)–(v) red:** the module moved. Park.
- **W2 Status 124:** see "Stop rule".
- No `PREDICTED RED` line.

### §9 item 231 — 2026-10-04 (chain H15a, re-pose of item 227, attempt 3 of 4)

**Outcome (2026-10-04 19:30 slot, implementer):** DONE — §10 H15a lands; `xl` entry 29 PENDING PREREQUISITE → READY. Module path-checked-out from `attempt/WF-7-step0e-r2-20261004T1714Z`, plus the one fix (`rss = _rss_summed_gib(comm)` on every rank before the rank-0 print); the binding grep found only `:348` (`:441` print-only, `:656` pure numpy/IO, `:711` all-rank `bcast`). W0 `20261005T003037Z_WF-7-H15a-r3-W0.log` Status 0, 14 passed, 60 s. W1 `20261005T003143Z_WF-7-H15a-r3-W1-orientation.log` Status 0, 131 s: 8 × `3 passed, 1 skipped`, exactly one `orientation control:` line (`:9795`) — 504 642 cells, φ̂ extent 2.816e-17 m (< 1e-12, asserted), **B = 128.39 s** (normal box; price gate ≤ 285 s passed). W2 `20261005T003402Z_WF-7-H15a-r3-W2.log` Status 0, 361 s, `[capture] rc=0`, one mesh line (binary green): (i) 507 266 cells (5.200e-3 vs band 0.01); (ii) reciprocity 1.041e-14; (iii) σ_max 0.999914648; (iv) 18 classes, worst spread 3.6997 %; (v) all 32 drives, max 3.328e-15 (stop rule did not fire: projection 466 s < 540 s); P17 reproduction 2.405e-11 / κ 5.996e-8; NEG 479.3×; npz 20 823 B, save + reload 0.014 s; **(L1)** all True; **(L2)** σ_max rel 2.596e-10, worst class-mean rel 4.279e-10 (band 1e-6); save control 3.147e-4 = 314.7× (inside [2e-4, 4e-4]); complex64 1.271e-8 (printed). Rule (i): module mtime 00:30:34Z < W1 start 00:31:43Z. Rule (a) as disclosed: the unset inline path was not run in-slot; its first run is entry 29.

*Drafted by a general-purpose agent for the 2026-10-04 18:00 daily review; the review reviews it.*

**`WF-7` step 0e, prerequisite half, re-posed (§10 H15a): item 227 exactly as written, plus one fix.** Item 227 died on a rank-0-only collective, not on price. The fix moves the summed-RSS `allreduce` out of the rank-0 print. The target is unchanged: on `FEM_EM_WF7_SAVE_NPZ`, the reloaded degree-1 32×32 reproduces σ_max and the 18 C16 × mirror class means within the imported `REUSE_REPRODUCTION_RTOL`; `xl` entry 29 goes READY.
- **Executor:** implementer. **Tier:** heavy, `-n 8`, complex (W1, W2); W0 standard, `-n 1`. **Class:** chain step (solve), measured median ≈ 23 min.
- **Cost:** ≈ 23 measured / ≈ 32 predicted min (15 fixed + ≈ 3 fix and grep + W0 ≈ 1 + W1 ≈ 2.5 + W2 ≈ 5.5 + ≈ 5 spare). Measured basis: W0 60 s and W1's build 120.57 s on 10-04 12:00 (below).
- **Depends on:** nothing. `git diff --stat b6f58e90 main -- tests/validation/test_wf7_f_human_ring_matrix.py` is empty (main's last touch is still `bfa8232c`), and the r2 branch's merge base is `b6f58e90`. Re-check this before starting.
- **Family:** `WF-7` step 0 (0, 0b, 0c, 0d, 0e). The last gate is 0d (`xxl` 2026-09-26). Since then: item 213 = 1 (price), item 227 = 2 (deadlock), **this item = attempt 3 of 4.** A failure here leaves one re-pose. **A fourth failure freezes the family** (daily-review.md, family cap): the review must then re-scope at a higher altitude as a new numbered step, or bank the negative. It may not reletter.

**Why item 227 died, measured** (`20261004T170232Z_WF-7-H15a-r2-W1-orientation.log`, Status 124 at 590 s).
- Ranks 1–7 each printed `3 passed, 1 skipped in 131.6 s` (`:9818–9836`).
- Rank 0 sat in `_rss_summed_gib(comm)`, an `allreduce`. The branch calls it inside `if comm.rank == 0:` in `test_wf7_default_orientation_control` (`:348–353` at `57789717`), where no other rank enters it.
- The default build completed: `mesh wall time 120.57 s` (`:9794`) ⇒ **B ≈ 121 s, a normal-speed box.**
- No orphans (`20261004T171254Z_WF-7-H15a-r2-cache-sweep.log`).
- W0 at `-n 1` cannot see this class of defect. W1 is the first window that can.

**Start from the branch.** `git checkout attempt/WF-7-step0e-r2-20261004T1714Z -- tests/validation/test_wf7_f_human_ring_matrix.py`. This is a module path checkout, never a record. `git diff --stat` must show that one file only. Everything item 213 and item 227 bound stays **as parked**: the save code, (L1), (L2), the controls, both knobs, both contract tests and the helper.

**The fix (tests only, one edit).** In `test_wf7_default_orientation_control`, compute the reduction on every rank before the branch:
- write `rss = _rss_summed_gib(comm)` straight after the `_default_orientation_extent(...)` call;
- inside `if comm.rank == 0:`, print `{rss:.3f}` in place of the call;
- nothing else changes, including the line's text and the assert.

**The grep (binding, before W0).** Run `grep -n "comm.rank == 0" tests/validation/test_wf7_f_human_ring_matrix.py` and read every block it hits for a `comm.` call or a helper that makes one. The drafting agent's reading at `57789717` found three hits:
- **`:348`:** the defect.
- **`:441`:** `say()`'s body only prints. Its f-string arguments are evaluated on every rank before the call, so `say(f"… {_rss_summed_gib(comm)} …")` at `:607–608` and `:809–810` is safe.
- **`:656`:** the save block. It calls `os.makedirs`, `np.savez` and `_reciprocity_ratio(s)`, which is pure numpy (`test_port_birdcage_ring_column.py:776`).

The `bcast` / `Barrier` at `:676–677` and `:711` run on all ranks. If the slot finds any other collective under a rank-0 branch, fix it the same way and list it in the commit.

**Anchors, negative controls, predictions, scope.** All are item 227's, unchanged: (i)–(v), the P17 reproduction, the column-scale control (asserted, ceiling ≈ 310×, 10-04 measured 314.7×), (L1), (L2) at `REUSE_REPRODUCTION_RTOL` 1e-6, the asserted mesh-level control (φ̂ extent < 1e-12 m), the complex64 control (printed), and W2's one-mesh-line binary. This item adds two checks.
- **Added binary (W1):** 8 `3 passed, 1 skipped` lines and **exactly one** `[WF-7 step1] orientation control:` line. A missing rank-0 line means the fix did not reach the code.
- **W1 predictions:** 504 642 default cells, extent ≈ 3e-17 m, build 120–133 s, wall 140–160 s.

**Windows.** Bash `timeout: 660000`, foreground. Every harness call is the bare, whole command, with `-s` everywhere.
- **W0**, standard, `-n 1`, predicted 14 passed, ≈ 60 s:
  `scripts/testing/run_and_log.sh WF-7-H15a-r3-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_save_knob_contract tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_orientation_knob_contract -v -s --tb=short'"`
- **W1**, heavy, `-n 8`. It runs the control alone and is the box-speed probe:
  `scripts/testing/run_and_log.sh WF-7-H15a-r3-W1-orientation "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_ORIENTATION_CONTROL=only timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short'"`
- **Price gate (after W1, binding, kept from item 227).** Read W1's printed `build_s`, B.
  - **B ≤ 285 s** ⇒ run W2.
  - **B > 285 s** ⇒ do not start W2. Park with W0 / W1 logs and the measured B: BLOCKED, "price — box > 2.2× slow".
- **W2**, heavy, `-n 8`, §5.1 durable capture with the trailing `; exit $rc` verbatim. Predicted 290–330 s:
  `scripts/testing/run_and_log.sh WF-7-H15a-r3-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf7 && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_ORIENTATION_CONTROL=split FEM_EM_WF7_SAVE_NPZ=/workspace/logs/wf7/wf7-s32-degree1-n8.npz timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short > /workspace/logs/wf7-h15a-r3-w2-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/wf7-h15a-r3-w2-raw.log; cat /workspace/logs/wf7-h15a-r3-w2-raw.log; exit \$rc'"`
- **Rule (i):** no module edit after W1 starts. Before committing, check the module's mtime against W1's start.
- **Rule (a):** as item 227. The unset (`inline`) path is not run in-slot, and its first run is entry 29's `xl` window. Disclosed; the review rules.
- **Stop rule:** item 227's. One W2 per slot, never a longer timeout, never `-n 12`, never a third heavy window. After any Status 124, run `pgrep -c python3` first, then the cache sweep only if the count is 0.

**Traps** (item 213's and 227's, plus):
- Any `comm.` call, or any helper that makes one (`_rss_summed_gib`, `_sheet_centre_z`, `_sheet_extent_along`), inside `if comm.rank == 0:` is a deadlock. `-n 1` cannot show it.
- Do not wrap a `say()` argument in a lambda or a rank-0 conditional: its collectives are safe only because every rank evaluates the f-string.

**Status it can move.**
- `docs/testing/xl-pending.md` entry 29: **PENDING PREREQUISITE → READY**, on W1 and W2 green. Edit the Status line only, in the same commit. Cite W2's log and the knob's value convention, plus the clause that `FEM_EM_WF7_ORIENTATION_CONTROL` unset is the landed inline path.
- The slot does not write, edit or queue entry 29's command, and does not touch `xl-queue.d/`.
- The §7 `WF-7` row stays ✅ and records the two knobs (records only).
- §10 H15a lands, so H15b becomes the review's to queue.

**Negative result.** In every case the slot parks on `attempt/`, marks BLOCKED in the same commit (rule (d)), and entry 29 stays PENDING.
- **W0 red:** fix once in-slot before W1. If it is still red, the slot parks; the review re-poses the contract.
- **W1 Status 124, or a missing rank-0 line:** the slot parks; the review reads the log for a second rank-0 collective against price.
- **W1 red** (the default sheets have φ̂ extent): the slot parks with a known-issues entry; the review rules on the generator. W2 does not run.
- **Price gate trips:** the slot parks; the review applies item 227's stop rule. That means W2's command as an `xl` priced-family variant of entry 21, or a quiet box, never a reletter.
- **(L1) red:** the slot parks; the review rules the save path.
- **(L2) red with (L1) green:** the slot parks and reports all 18 deltas; the review routes re-recording to `record-reconciler` under (1*), never in-slot.
- **Any (i)–(v) red:** the slot parks; the review rules on what moved.
- Any of these spends attempt 3; one remains before the freeze.
- No `PREDICTED RED` line.
