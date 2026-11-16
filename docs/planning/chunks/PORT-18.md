# PORT-18 — Does the lumped-sheet port read its current off one side of a discontinuous component

 **Does the lumped-sheet port read its current off one side of a discontinuous component — and is that side assigned differently on the four C4 sheets?** — opened 2026-09-10 03:00 review on `GEO-31`'s verdict (same support, same boundary, **different interior cut** at every rung, including the ×1 rung every port gate passed on). **The mechanism this row tests, read off the code rather than inferred:** on the four-port birdcage the sheet `SHEET_IFACE + i` is the interface between `PORT_LOWER + i` and `PORT_UPPER + i` — the gap box's mid-plane, so its normal is (to be verified in-probe) `±ẑ` — and the fixture drives it with `drive_direction=(0, 0, 1)` (`tests/validation/test_port_birdcage_four_port.py:397`), i.e. **`ĥ ∥ n̂`**. `sheet_terminal_current` reads `∫ E('+')·ĥ dS` and `lumped_port_linear_term` injects `K_imp·v('+')` (`src/fem_em_solver/ports/lumped.py:306, 333`); `_restrict`'s docstring justifies `'+'` because `n̂ × E` is single-valued (Jin (1.60)) — true of the **tangential** trace the bilinear term uses, **not** of the **normal** component `E·n̂`, which N1curl leaves discontinuous across a facet. Which cell DOLFINx calls `'+'` is a per-facet property of the cut, so four congruent patches cut differently can read four different mixes of lower-side and upper-side traces. That would put the asymmetry into both the source and the read-back on every rung, amplified as the cut coarsens relative to the field — the frame `GEO-31` reading (i) set. Executor **`mesh-probe`**, measurement-only, 🧪 by the §3 rule, **no solve**, real build. **Import `scripts/probes/geo30_birdcage_c4_refinement_census.py`'s `RUNGS` / `_build_rung` and `GEO-31`'s `SHEET_IFACE + i` interface-tag selection verbatim** (as `GEO-31` did); new `scripts/probes/port18_sheet_plus_side_census.py`. **The measurement**, per rung (all four `RUNGS`), per sheet: (0) `∫ (n('+')·ẑ)² dS / A_i` — the premise check, printed; (1) the **`'+'`-side lower-half fraction** `f_i = ∫ χ_L('+') dS / A_i` with `χ_L` the DG0 indicator of cell tag `PORT_LOWER + i`; (2) for a smooth, C4-invariant field **not** in N1curl degree 1, `E = ẑ·cos(2π(z − z̄)/λ)·(1 + (x² + y²)/ρ₀²)` (`z̄` the common sheet-centroid height, `λ` = 4× the sheet's gap height, `ρ₀` the leg radius, all measured in-probe and printed) interpolated into N1curl degree 1: the read-backs `R_i^+ = ∫E('+')·ẑ dS / A_i`, `R_i^−` and `R_i^avg` at pinned quadrature degrees `q ∈ {2, 4, 8}`, the expression comparand `∫ g dS / A_i` at `q = 8`, and the four-sheet spread `(max − min)/|mean|` of each. (3) Repeat (1) at `-n 1` on the ×1 and ×0.75 rungs only — side assignment can depend on the partition. **Anchors (asserted — closed forms and an identity):** (A) partition identity `∫ (χ_L + χ_U)('+') dS = A_i` to ≤ 1e-12 relative on every sheet at every rung — without it (1) measures a broken indicator; (B) interpolation-exactness: `E_lin = (−z/ρ₀, 0, 1 + x/ρ₀)` is of the form `a + b × x`, so N1curl degree 1 represents it exactly and its trace is continuous: assert `R_i^+ = R_i^− = 1 + x̄_i/ρ₀` to ≤ 1e-9 relative on every sheet at `q = 2`, with `x̄_i = ∫x dS / A_i` measured in-probe. **Negative control (asserted — arithmetic, not predicted):** the same `E_lin` is deliberately *not* C4-invariant; for four sheets 90° apart at radius `ρ₀`, `Σ cos θ_i = 0` and `max − min` of `cos θ_i` is ≥ √2, so its four-sheet spread must read **≥ 1.40** — a spread metric that returns ≈ 0 regardless fails here. **Printed, asserted nowhere (rule (e)):** the `f_i` table at both rank counts; the spreads of `R^+`, `R^−`, `R^avg` and the expression by `q` at all four rungs side by side with the imported gate spreads (`ANS-4` ×1 `Z` classes 0.1012 / 0.0916 / 0.0654 %, ×0.75 0.5390 / 0.4591 / 1.6886 %, known-issues 🔴 2026-09-09); and the discriminator sentence: **(a)** `f_i` equal across the four sheets at every rung and `R^+` spread ≲ the `R^avg` spread ⇒ the read-back's side choice is excluded, and the solve (the source term's `v('+')·ẑ` included) is the remaining suspect; **(b)** `f_i` differs across sheets and the `R^+` spread exceeds the `R^avg` spread and tracks the rungs that broke ⇒ the `'+'` restriction of a normal component is a port-model defect in `lumped.py`; **(c)** spreads fall with `q` at fixed `f_i` ⇒ quadrature, curable. **Tier / ranks / cost:** four meshes (`GEO-31` measured 144 s for exactly these at `-n 2`, meshing 25–37 s each) plus DG0/N1curl interpolation and a few dozen facet assemblies ⇒ window 1 `-n 2`, `timeout -k 30 500`; window 2 `-n 1` on two rungs, `timeout -k 30 300`; `-s` on both, real build; standard by expectation, heavy by ceiling. **Traps already paid for:** every `SpatialCoordinate`-bearing facet form pins `metadata={"quadrature_degree": q}` (`POST-5` step 1, nine-minute FFCx compiles, poisoned cache entries); interior-facet forms need a restriction on every argument (`_restrict`); no UFL ordering comparison and no `sqrt` inside one (`OPS-22`, `WF-6` step 4c) — the field above uses neither; the DG0 indicator must be set on **ghost** cells too, or `'+'` on a partition-boundary facet reads 0 (set owned values, then `x.scatter_forward()`); `cell_tags.values`, `assemble_scalar` and centroids are rank-local — `allreduce` before printing; `timeout -k 30`; no `run_in_background`; `-k a or b` splits inside a quoted container command; sweep `find /root/.cache/fenics -name '*.c' -size 0` and delete stubs only. **Scope:** measures the side assignment and the read-back only, not the solve; closes nothing and rehabilitates no rung. `ANS-4` keeps ✅ + INCONCLUSIVE, `WF-6` stays 🟡, `ANS-4` step 2a / `WF-6` step 4f stay blocked, and no port band, record or `src/` line moves. **Negative result:** if (0) reads far from 1 the premise `ĥ ∥ n̂` is wrong: print it and stop. If (b), open a known-issues row naming `sheet_terminal_current` and `lumped_port_linear_term` and stop. **Do not switch to `ufl.avg` in-slot:** it would move every `PORT-9` / `PORT-11` / `PORT-16` / `ANS-4` / `WF-6` record at once, and that is a review's ruling with its own re-run plan. If (a), record it and stop: the finding moves the suspect onto the solve. **MEASURED 2026-09-10 04:30 slot — PREMISE (0) FALSE; stopped on the row's first negative branch** (`mesh-probe`, `scripts/probes/port18_sheet_plus_side_census.py`; window 1 `-n 2` all four rungs Status 0 / 154 s, `20260910T093526Z_PORT-18.log`; window 2 `-n 1` ×1 + ×0.75 Status 0 / 68 s, `20260910T093810Z_PORT-18.log`). **Asserted self-tests green in both windows:** (A) partition identity worst 1.089e-15 vs 1e-12; (B) `E_lin` read-backs `R^±` = `1 + x̄_i/ρ₀` to ≤ 9.99e-16 (relative denominator `max(|1 + x̄_i/ρ₀|, 1)`, since the θ = 180° sheet's target is ≈ 0 — disclosed in the docstring); negative control spread **2.000000** vs ≥ 1.40 (`…093526Z…:7011–7025, 7144`). **Premise (0):** `∫(n('+')·ẑ)² dS / A_i` = **0.000000000000000 on all 16 sheet readings at `-n 2` and all 8 at `-n 1`** (`:7008, 7037, 7066, 7095`; `…093810Z…:3517, 3546`), corroborated independently of the form by the sheet bounding boxes — extents (1.400e-02, ≈1e-18, 8.000e-03) m on P1/P3 and (≈1e-18, 1.400e-02, 8.000e-03) m on P2/P4 (`:7007`): **each sheet is a vertical radial–axial plane with an azimuthal normal**, not the gap box's horizontal mid-plane this row read off the code. So `ĥ = ẑ` lies **in** the sheet, `E('+')·ẑ` is a **tangential** component N1curl keeps single-valued, and the mechanism this row tests (a cut-dependent `'+'` side on a discontinuous normal component) **does not exist on this fixture** — for the read-back and, by the same argument, for the source term's `v('+')·ẑ`. The (a)/(b)/(c) predicates were not evaluated; no known-issues row, no `ufl.avg`, no `src/` or band change. **Printed, asserted nowhere, and consistent with the premise reading:** `f_i` (here the `PORT_LOWER + i` box half, which on this geometry is not a z-half) varies across sheets *and* rearranges completely with the partition (×1: 0.714/0.148/0.091/0.722 at `-n 2`, 0.544/0.959/0.017/0.193 at `-n 1`), yet `R^+`, `R^−`, `R^avg` are **identical to every printed digit** at `q` = 2, 4 and 8 and at both rank counts — side assignment and quadrature both demonstrably inert. **One lead for the review, printed only:** the smooth field's N1curl-interpolant read-back spread tracks the rungs that broke — ×1 7.223213e-06, ×0.75 3.520702e-04 (≈ 49×), `resolution` 0.012 3.139712e-06, 0.0095 7.307614e-04 (≈ 233×) — while the exact expression comparand sits at ≤ 4.2e-15, and the sheet pairings match `GEO-31`'s centroid readings (×1 P1 = P2, P3 = P4; ×0.75 P1 = P3, P2 = P4; 0.012 P4 alone; 0.0095 all four). The magnitudes are 0.035% / 0.073%, below the imported `ANS-4` ×0.75 `Z` class spreads (0.5390 / 0.4591 / 1.6886 %), so the read-back of a *given* field is cut-sensitive but not by itself at gate scale; whether the solve amplifies it is not measured here. Disposition is the review's. **Disposition, 2026-09-10 10:30 review:** accepted as measured (`20260910T093526Z_PORT-18.log:7008, 7037, 7066, 7095`; `20260910T093810Z_PORT-18.log:3517, 3546`; both asserted self-tests green; 🧪, no audit owed). **The wrong premise was the 03:00 review's code reading, not the executor's.** The generator's own comment says the sheet plane contains the leg axis *and* the radial direction, so its normal is azimuthal (`src/fem_em_solver/io/mesh.py:3972–3974`). Premise check (0) caught the error in 222 s of compute. The `'+'`-side mechanism is **excluded** for both the read-back and the source term. The executor's relative denominator in (B), `max(|target|, 1)`, is ratified: the θ = 180° sheet's target `1 + x̄/ρ₀` is ≈ 0, and the absolute residual ≤ 9.99e-16 is the stronger statement either way. The printed lead is taken up by **`GEO-32`**: the interpolant read-back spread tracks the broken rungs with `GEO-31`'s sheet pairings, and a C4-congruent cut is the one change that removes that lead without a solve. `ANS-4` step 2a and `WF-6` step 4f stay blocked. | 🧪 **MEASURED 2026-09-10 04:30 slot** — premise false; measurement-only, `🧪` by the §3 rule

### §9 item 246 — 2026-10-05 (lane 1)

*Drafted by a general-purpose agent for the 2026-10-05 03:00 daily review; the review reviews it.*

**`PORT-18` step 2: the sheet's two one-sided current readings on the exact rectangles.** This is the lane-1 target the 10-03 weekly named (PROJECT_PLAN.md:4643–4644): "on the exact-rectangle fixture the sheet's two one-sided current readings agree at the imported band ⇒ ✅".
- **Executor: implementer, not `mesh-probe`.** The item asserts, and it runs one solve.
- **Tier:** standard by prediction, heavy by ceiling (`timeout -k 30 300`). Complex build. Two windows: W1 at `-n 2`, W2 at `-n 1`. Side assignment rearranges with the partition (09-10: `f_i` at ×1 read 0.714/0.148/0.091/0.722 at `-n 2` and 0.544/0.959/0.017/0.193 at `-n 1`), so both widths are needed.
- **Class:** lane 1. The class median is ≈ 8 min (restock 10-04 18:00). **Cost:** ≈ 8 measured / ≈ 35 predicted min (15 fixed, ≈ 12 code, two windows of ≈ 2–3 min each, spare).
- **Depends on nothing.** No `src/` change: one new test module only.

**What "the two one-sided current readings" are, read off the code.**
- `sheet_terminal_current` (`src/fem_em_solver/ports/lumped.py:338–364`) returns `I⁺ = (1/(R h)) ∫_S (E('+') + E_src)·ĥ dS`. It restricts through `_restrict` (`:270–281`), which hard-codes `'+'` on an interior sheet.
- The other reading is `I⁻`: the same integrand with `'−'`. The 09-10 probe measured exactly this pair as `R⁺` / `R⁻`. For the smooth field see `scripts/probes/port18_sheet_plus_side_census.py:258–259`; for `E_lin` see `:283–284`.
- 09-10 measured the sheet normal as azimuthal and `ĥ = ẑ` as lying in the sheet (`20260910T093526Z_PORT-18.log.gz:7007–7008`). So `E·ĥ` is a **tangential** component. N1curl keeps that single-valued, which makes `I⁺ = I⁻` an identity up to round-off whichever cell DOLFINx calls `'+'`. The item asserts that identity on the sheet the ports actually use, for a solved field.

**The imported band (never restated).**
- `EXACTNESS_BAND` = 1e-9 (`scripts/probes/port18_sheet_plus_side_census.py:84`). This is the band that probe held `R^±` to in (B) (`:493–497`).
- `PREMISE_BAND` = 1e-6 (`:86`), the probe's own threshold on the premise reading `∫(n('+')·ẑ)² dS / A`.
- Import both with the probe's own path-insert idiom (`:71–77`). Precedent for a test importing from `scripts/probes`: `tests/validation/test_post8_owner_rule_ans2.py`.
- If the import fails or has side effects, **stop BLOCKED**. Do not retype the numbers.

**The comparison from the 09-10 record, and a correction to the target's wording.** The weekly's phrasing expects the exact sheet to collapse a ragged-sheet discrepancy. **The record holds no such discrepancy.**
- 09-10 read the **full gap face `F_i`** (58 triangles at ×1, `…093526Z_PORT-18.log.gz:7000`). `F_i` is a CAD face on both fixtures. The probe did not read the narrowed strip `N_i` that the ports integrate over, and it read no solved field.
- On `F_i`, `R⁺ = R⁻` to every printed digit at `q` = 2 / 4 / 8 (13 significant digits, i.e. ≤ ≈ 2.7e-13 relative, `:7012–7020`). `E_lin`'s `R^±` matched the closed form to ≤ 8.88e-16 (`:7023–7024`).
- On the old fixture the ragged part was `N_i`: +4.20e-2 in area (`20261001T215942Z_OPS-61-F10a-W1.log:3582`). Its one-sided readings were never measured.
- **This item therefore tests no collapse.** It carries the identity onto the exact `N_i` (40 triangles at both widths, `20261001T220051Z_OPS-61-F10a-W2.log:3591–3594`) and onto a solved field. The predicted reading is round-off (≈ 1e-15 … 1e-13), not a fall.
- **The review should say so in the row amendment below.**

**Code (tests only): new `tests/validation/test_port18_exact_strip_one_sided.py`, `complex_only`.**
- **Fixture (module scope).** Call `rung = _four_port_rung("PORT-18 exact 10 MHz", np.zeros(4), degree=1, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` (`tests/validation/test_port_birdcage_leg_offset_sweep.py:199`).
  - This is item 192's base call (`test_ans4_exact_strip_records.py:435–448`) at the default `FREQUENCY_HZ` = 10 MHz, with κ off (`width_correction_kappa=None`, the 10-02 ruling at PROJECT_PLAN.md:3692).
  - It hoists `create_connectivity` / `create_entity_permutations` (`:294–295`).
  - Then `solved = _solve_driven(rung, "P1")` (`tests/validation/test_birdcage_b1_plus_map.py:142–203`). It keeps `E` and returns the production `I⁺` for all four sheets in `solved["currents"]`.
- **`I⁻`.** Rebuild the sheets exactly as `_solve_driven` does (`[s.sheet(driven=(s.port_id == "P1")) for s in rung["specs"]]`). Mirror `lumped.py:356–364` with `r = lambda f: f("-")`, using `lumped._sheet_measure` for the measure. Allreduce, then divide by `h`.
  - Also build a **`'+'` copy** of the same mirror, for the fidelity check (B0).
- **Jump metric.** `J_i[g] = sqrt( ∫_{N_i} |(g('+') − g('−'))·ĥ|² dS / ∫_{N_i} |g('+')·ĥ|² dS )`.
  - Write the integrand as `ufl.inner(d, d)`, with `d` the scalar difference.
  - Assemble on every rank and allreduce. Assert `|imag| ≤ 1e-12·real` silently, print the imaginary part, and take `sqrt` in numpy.
  - No `sqrt` and no comparison inside UFL.
- **Interpolated fields** (the probe's definitions, `:243–254`, `:276–281`):
  - `E_lin = (−z/ρ₀, 0, 1 + x/ρ₀)` and the smooth field, both in N1curl degree 1 on the rung's mesh.
  - `ρ₀` is the mean `N_i` centroid radius, printed beside `RING_RADIUS`. `z̄` is the mean centroid height, and `λ` = 4 h.
  - Every `SpatialCoordinate`-bearing facet form pins `metadata={"quadrature_degree": q}`.
- **Planted control field.** `g = s·ẑ`, with `s = χ_L − χ_U` the DG0 indicator difference of `PORT_LOWER + i` / `PORT_UPPER + i`. Set owned values from `cell_tags`, then `x.scatter_forward()` (the probe's `_indicator` idiom, `:141–150`).

**Anchors (asserted, every sheet, both widths unless stated).**
- **(P0) Fixture identity.**
  - Cells vs `ANS4_EXACT_STRIP_CELL_RECORD` within `CELL_COUNT_BAND`.
  - Each `N_i` area vs `SPEC_AREA_M2` (`test_ans4_exact_strip_records.py:131–133`) within `STRIP_BAND` (`test_ops61_exact_strip_ports.py:58`).
  - **At `-n 2` only** (`ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`): `‖S − ANS4_EXACT_STRIP_S_P1_RECORDS["10"]‖/‖rec‖ ≤ STEP3A_RECORD_RTOL`. At `-n 1` it is printed.
  - Backed by item 192 W1 (`20261003T093408Z_ANS-4-item192-W1-p1.log:1899, 2117–2132`).
- **(P1) Premise.** `∫(n('+')·ĥ)² dS / A_i ≤ PREMISE_BAND` on every `N_i`, i.e. ĥ lies in the sheet.
  - *Asserted*: 24 of 24 readings on `F_i` were 0.000000000000000 (`…093526Z…:7008, 7037, 7066, 7095`; `…093810Z…:3517, 3546`). `N_i ⊂ F_i` is coplanar.
- **(A) `E_lin` closed form.** `R_i^± = 1 + x̄_i/ρ₀` within `EXACTNESS_BAND`, with the relative denominator `max(|target|, 1)` (ratified by the 09-10 10:30 review).
  - *Asserted*: ≤ 8.88e-16 on `F_i` (`:7023–7024`).
- **(B0) Mirror fidelity.** The test's `'+'` copy reproduces `solved["currents"][P_i]` to ≤ `EXACTNESS_BAND` relative. Without this, (B) compares a production reading with a different integrand.
- **(B) The weekly's target.** For the P1-driven solve: `|I_i⁺ − I_i⁻| / |I_i⁺| ≤ EXACTNESS_BAND` on all four sheets, and `J_i[E] ≤ EXACTNESS_BAND`. The same `J_i` bound holds for the smooth interpolant.
  - *Predicted* ≈ 1e-15 … 1e-13 (the 09-10 interpolant reading on `F_i`); printed beside the measured value.
  - The undriven sheets carry O(0.5)-of-driven currents (|S₂₁| ≈ 0.47 in the record), so no denominator is near zero. Print the four `|I_i⁺|`.

**Negative control (asserted: arithmetic, not predicted).**
- `J_i[g] = 2` exactly for the planted DG0 field.
  - Every `N_i` facet separates a `PORT_LOWER + i` cell from a `PORT_UPPER + i` cell. `F_i` is built from that cell-tag pair (`io/mesh.py:290–373`, via `_sheet_areas`), and `N_i ⊂ F_i` (`OPS-61` (A3)).
  - So `g('+') = −g('−') = ±ẑ`, whichever side is `'+'`.
- Assert `|J_i[g] − 2| ≤ EXACTNESS_BAND`. The separation from (B)'s band is 2e9. 2 is this field's closed form, not a tunable factor.
- A jump metric that is blind to a discontinuity (one that reads only `'+'`, for example) returns 0 here and fails.
- **Printed beside it:** the signed one-sided means of `g`, which are `±(2f_i − 1)` with `f_i = ∫χ_L('+') dS / A_i`, at each width. This shows the side assignment the identity is indifferent to.

**Printed, never asserted.**
- `I_i^±` to 10 digits, `|I⁺ − I⁻|/|I⁺|`, and `J_i` for each field.
- `f_i`.
- The smooth-field four-sheet spread on `N_i`, beside 7.223213e-06 on `F_i` ×1 (`…093526Z…:7012`). That is `GEO-32`'s lead, not this item's.
- Rung, sweep and solve times.

**Windows** (foreground, Bash `timeout: 660000`; each harness call is the whole Bash command, bare and repo-relative; `-s`; no pipes; no `-k` selection).

`scripts/testing/run_and_log.sh PORT-18-step2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_port18_exact_strip_one_sided.py -s -v --tb=short'"`

`scripts/testing/run_and_log.sh PORT-18-step2-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_port18_exact_strip_one_sided.py -s -v --tb=short'"`

**Cost, from logs.**
- Exact ×1 rung: 29.77 s, and four degree-1 driven solves in 9.00 s, at `-n 2` (`20261003T093408Z_ANS-4-item192-W1-p1.log:1899`). The rung is 27.05 s at `-n 1` (`20261001T215942Z_OPS-61-F10a-W1.log:3575`).
- `tests/environment` plus light port modules: 51 s at `-n 1` (`20261002T201726Z_PORT-3-W2.log`).
- Predicted W1 ≈ 100–150 s and W2 ≈ 110–180 s, inside the 300 s ceiling.
- **Box-speed gate:** if the rung print exceeds 60 s at `-n 2`, the box is slow. That is a price, not a defect. A window killed at 300 s stops the item BLOCKED with the reading; never lengthen the timeout.

**Traps.**
- **Rank-local data:** `cell_tags.values` and `facet_tags.values`. The DG0 indicator must reach ghost cells (`scatter_forward`), or `'+'` on a partition-boundary facet reads 0. `assemble_scalar` is reduced before any assert.
- **Collectives:** every rank runs every collective, and asserts come after them. A `comm.` call inside `if comm.rank == 0:` deadlocks (item 227).
- **Complex build:** `ufl.inner` conjugates its second argument, so (B)'s difference is taken on the restricted values before `inner`. No UFL ordering comparison and no `sqrt` (`OPS-22`, `WF-6` 4c). Take real parts only after the reduction, and print the imaginary parts.
- **Compiles:** pin `quadrature_degree` on every `SpatialCoordinate` facet form (`POST-5`). Sweep orphan `.c` stems only with `pgrep -c python3` reading 0, as the rubric's one-line sweep. After a kill, run `pgrep -c python3` first.
- **The S premise** is asserted at rank width 2 only. At `-n 1`, print it.
- **Run hygiene:** rule (i); no edit after the last window.

**Scope.**
- Degree 1, 10 MHz, P1 drive, on the exact ×1 mesh only. Degree 2, the production default, is **not** exercised. The identity is degree-independent in theory but unmeasured there.
- The source term's `v('+')·ĥ` rests on the same tangential argument and is not measured.
- No ragged rung is re-tested or rehabilitated. `ANS-4` step 2a and `WF-6` step 4f stay as they are.
- No band, record, κ setting or `src/` line moves.
- A ✅ licenses **no** `S₁₁` / `Z_in` / `f₀` figure. It says only that the read-back's `'+'` restriction is not a source of asymmetry on exact-rectangle sheets.

**Negative result.**
- **(B) red** (with (P0), (P1), (A), (B0) and the control green): the tangential trace disagrees across a facet. That is an orientation or permutation defect, not physics.
  - Open a known-issues row naming `sheet_terminal_current` and `_restrict`.
  - Park on `attempt/PORT-18-step2-<UTC>` and mark BLOCKED in the same commit (rule (d)).
  - **No `ufl.avg` switch in-slot.** It would move every port record at once; this is the row's own instruction.
- **(P1) red:** the exact-strip sheet is not the plane the 09-10 record measured. Stop and report.
- **(P0) red:** the fixture drifted. On an image change this is `record-reconciler` work; otherwise stop.
- **(B0) or control red:** the instrument is wrong. Fix it once in-slot, disclose, and re-run the window. If it is still red, BLOCKED.
- There is no PREDICTED RED line; the drafter predicts green.

**Status it can move.** `PORT-18` 🧪 → ✅.
- **Row amendment (rule (j)), in the review's own commit.** The §7 row (PROJECT_PLAN.md:2683) carries no Done-when, only the 09-10 measurement. The review appends:
  > **Done-when (2026-10-05 03:00 review, rule (j); the 10-03 weekly's lane-1 target, §9 item 246):** on the exact-rectangle ×1 fixture (`leg_sheet_strip_fraction = GATED_WIDTH_FRACTION`, κ off), for a P1-driven degree-1 10 MHz solve at `-n 2` and `-n 1`, `sheet_terminal_current`'s `'+'` reading and its `'−'` mirror agree within `EXACTNESS_BAND` (probe:84) on all four sheets, together with the jump metric `J_i[E]`. `E_lin` meets its closed form, the premise ĥ ⊥ n̂ holds, and the planted DG0 discontinuity reads `J = 2` ⇒ 🧪 → ✅. *(The 09-10 record holds no ragged-sheet one-sided discrepancy to collapse: it read the CAD face `F_i`, where `R⁺ = R⁻` to print precision; this closes the identity on the gated strip `N_i` and on a solved field.)*
- The Tier cell also gains: "step 2: standard by prediction, heavy by ceiling, `-n 2` + `-n 1`, one degree-1 solve".

**§9 line 2:** → docs/planning/chunks/PORT-18.md#§9-item-246 

**Outcome — ✅ DONE 2026-10-05 09:00 CDT slot (implementer; `PORT-18` 🧪 → ✅).** New `tests/validation/test_port18_exact_strip_one_sided.py`, no `src/` change; `EXACTNESS_BAND` / `PREMISE_BAND` imported from the probe by its own path-insert idiom (no side effects). W1 `-n 2` 18 passed, Status 0, 121 s (`20261005T140253Z_PORT-18-step2-W1.log`); W2 `-n 1` 18 passed, Status 0, 117 s (`20261005T140501Z_PORT-18-step2-W2.log`). Rung 34.50 s at `-n 2` (box-speed gate 60 s: normal box), 41.00 s at `-n 1`.

**Step 2 narrative (measured 2026-10-05).**
- **(P0)** cells 117 206 = record at both widths; every `N_i` 40 facets, area vs SPEC ≤ 1.6e-15; `‖S − S_rec(10)‖/‖S_rec‖` = 6.601e-11 at `-n 2` (asserted, band 1e-6), 6.598e-11 at `-n 1` (printed) (`…W1.log:1919–1922`, `…W2.log:1897–1898`).
- **(P1)** `∫(n('+')·ĥ)² dS / A_i` ≤ 1.2e-33 on all 8 readings — ĥ = ẑ lies in the exact strip, as on `F_i` in step 1.
- **(A)** `E_lin` `R^±` vs `1 + x̄_i/ρ₀` ≤ 1.78e-15 (ρ₀ = 7.000000000000e-02 m = `RING_RADIUS`, z̄ ≈ 1e-20 m, h = 8.0e-3 m, λ = 3.2e-2 m).
- **(B0)** the `'+'` mirror reproduces the production `I⁺` to 0.000e+00 on all four sheets at both widths.
- **(B)** `|I⁺ − I⁻|/|I⁺|` ≤ 2.1e-16 at `-n 2` and ≤ 5.0e-16 at `-n 1` (predicted 1e-15 … 1e-13, so at or below the low end); `J[E]` ≤ 3.25e-15; `J[smooth]` ≤ 3.37e-15. `|I⁺|` = 1.3918e-2 (P1, driven), 4.7362e-3 / 4.4678e-3 / 4.7360e-3 A (P2–P4); no denominator near zero. Numerator imaginary parts read exactly 0.
- **Negative control** `J[g]` = 2.000000000000000 on every sheet at both widths (|J − 2| = 0.000e+00).
- **Printed only.** `f_i` rearranges with the partition, as in step 1: 0.326 / 0.040 / 0.083 / 0.029 at `-n 2`, 0.717 / 1.000 / 0.035 / 0.227 at `-n 1`, and the signed one-sided means of `g` read ±(2f_i − 1) as predicted. The identity is indifferent to that side assignment. The smooth interpolant's four-sheet spread on `N_i` reads 1.96e-15 (`-n 2`) / 2.20e-15 (`-n 1`) against 7.223213e-06 on `F_i` ×1 in step 1. On the exact strips the cut-sensitivity lead `GEO-32` took up is at round-off. This is printed, not asserted, and this item draws no conclusion from it.
- **Scope held.** Degree 1, 10 MHz, P1 drive, exact ×1 only. Degree 2 and the source term's `v('+')·ĥ` are not measured. No band, record, κ or `src/` line moved.
