# OPS-84 — the (S2-tip) comparand's form on a field-derived source: I41's matched rod pair

*History file for the §7 `OPS-84` row (`OPS-46`). Opened by the 2026-11-06 10:30 daily review (lane 1, known-issues diagnosis) from the known-issues entry of 2026-11-06, whose Cause line `OPS-83` filled on the stub and left one clause open: "**Not measured:** that the clinical-scale reds (entry 55, U127) are wholly this form gap; the stub shows the mechanism at 1e-3, not its size on a degree-1 \|E\|²." The weekly re-poses the gate; this chunk only measures.*

**Why.** `OPS-83` (`20261106T120232Z_OPS-83-W1.log:78–84`) showed the mechanism on a hand-made checkerboard: R_DG0 is Δt-independent at 1.6e-3 / 2.1e-3, and against the consistent-mass projection R_Π falls to 7e-6 and scales with Δt. The reds that matter are on field-derived sources — the DG0 cell average of ½σ\|E\|² from a degree-1 N1curl solve near a PEC rod tip, where the source is far rougher than a factor-3 checkerboard. U127 reads 4.3e-3 / 5.9e-4 (P1 drive, exact-rectangle rod mesh) and entry 55 reads 2.5e-3 / 6.4e-3 (F-clinical), with the gel twin at 3.1e-3 / 1.3e-3. Review ruling 5 on U127 (`WF-15.md`, 2026-11-05 03:00) also named a conduction term on a tip ball, ≈ 1e-5 to 1e-3 in review arithmetic, which is O(Δt) and therefore shows in R_Π, not in its Δt-independent part. Neither U127's mesh (its `Study` route lives only on `attempt/WF-15-s14-20261105T212843Z`) nor entry 55's (`xl`) can be read in a slot. I41's matched pair can: it is on `main`, it carries the rod and its gel twin on one mesh family with the same tip balls (`TIPS`, 1 g), and its last window ran 310 s at `-n 2` (`20261105T183330Z_OPS-81-W2-i41.log:3772`, `:3798`).

**Done-when.** `tests/validation/test_ops84_s2tip_form_rod_pair.py` lands, every asserted gate below green at `-n 2` in the complex build, the separating table and the pre-stated reading printed, and the known-issues entry's Cause line gains the field-derived reading (the "Not measured" clause answered, confirmed or refuted); the entry stays open and resolves with the weekly's re-pose ⇒ ✅.

### §9 item 860 — 2026-11-06 (lane 1, known-issues diagnosis; implementer; heavy, `-n 2`)

**Code.** New `tests/validation/test_ops84_s2tip_form_rod_pair.py`, complex build (`@complex_only`), knob `FEM_EM_OPS84=1` (module fixture skips without it, as I41's does at `test_imp3_rod_heating_matched.py:202–206`). No `src/` edit; no landed module edited.
- **Imports, unchanged:**
  - from `tests.validation.test_imp3_rod_heating_matched` (I41): `_arm_with_403`, `_merge_rod_into_phantom`, `_census`, `ROD_CELLS`;
  - from `tests.validation.test_imp3_rod_heating_phantom` (I39): `_outer_robin_tags`, `_parent_wall_count`;
  - from `tests.validation.test_imp3_rod_birdcage`: `_arm` (as `_rb._arm`, the way I41 calls it), `TIPS`;
  - from `tests.mesh.test_birdcage_implant_rod`: `_build_rod`, `_build_kept_rod`, `KEPT_ROD_CELL_RECORD`;
  - from `tests.validation.test_mat15_coil_driven_heating_phantom`: `_transient`, `_tag_submesh`, `_cell_average_dg0_tag`, `DT_1`, `RHO_C`, `T_INF`, `S2_BAND`, `E_BAND`, `_max_imag`;
  - from `tests.validation.test_mat10_sar_to_heat_sphere`: `_to_sub`, `_allsum`, `B_BAND`, `ROBIN_TAG`;
  - from `tests.validation.test_mass_averaged_sar_standard_masses`: `ONE_GRAM_KG`, `QUADRATURE_DEGREE`;
  - from `tests.validation.test_birdcage_sar_mass_averaged`: `PHANTOM_RHO_KG_PER_M3`;
  - from `tests.validation.test_port_birdcage_lumped_column`: `PHANTOM_CELL_TAG`;
  - from `tests.validation.test_ops83_s2tip_form_on_stub`: `_projections`, `CONST_RATIO`, `SCALING_RATIO`;
  - `fem_em_solver.post.sar`: `mean_sar`, `averaging_ball_radius`.
- **Read every signature first** (`_arm` at `test_imp3_rod_birdcage.py:295`; its hook receives `(label, rod, msh, cell_tags, e_of, sigma_field, rho_field)`, and `e_of` holds `"P1"` … and `"ccw"` / `"cw"`, `:340–341`).
- **Run order = I41's `_run` without its control arm** (`test_imp3_rod_heating_matched.py:133–178`): build R (`_build_rod`) and M (`_build_kept_rod`), assert (A0m) exactly as I41 does (`:149–151`); arm M first through `_arm_with_403(built_m, …, map_gel=True)` with a hook on the merged tags (`_merge_rod_into_phantom`); then arm R through `_rb._arm(built_rod, …, rod=True)`. **All thermal work runs inside the hook** — the fields are freed when `_arm` returns.
- **In each hook, a new `_slope_arm(label, msh, cell_tags, e_of, sigma_field, G_ccw)`**, for each source s ∈ {`"ccw"`, `"P1"`}:
  1. P_s = `mean_sar(e_of[s], …, subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]` (real part), as I39's `_heat_arm` (`test_imp3_rod_heating_phantom.py:153–155`). G: on ccw, arm M's G = 1/P_ccw(M) and arm R reuses it (I41's scheme, so the reproduction control below is like for like); on P1, G = 1/P_P1 per arm (the slopes are scale-free).
  2. `sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)`, `ft, counts = _outer_robin_tags(sub)` once per arm; q_s = `_to_sub(*_cell_average_dg0_tag(msh, cell_tags, PHANTOM_CELL_TAG, G·½σ·inner(e, e)), sub, parent_ids)` exactly as `:168–171`; (B) ∫_sub q_s = G·P_s.
  3. BE one step from T∞: `_transient(sub, ft, q_s, dt=Δt, n_steps=1, theta=1.0)` for Δt ∈ {`DT_1`, `DT_1`/10} — four solves per arm, eight in all. Global (S2) per `:188–190` with `G·P_s / V_sub`.
  4. `pq, plq = _projections(q_s, V1, tag)` with V1 = `res.temperature.function_space` and **a tag unique per arm and source** (`f"ops84_{arm}_{s}"`; it is the PETSc options prefix).
  5. **Ball reads** at both `TIPS` with a = `averaging_ball_radius(ONE_GRAM_KG, PHANTOM_RHO_KG_PER_M3)`, through I39's indicator form (`:202–216`: `ufl.conditional(ufl.lt(ufl.real(…), ufl.real(rad_c) ** 2), 1, 0)`, `dx` with `metadata={"quadrature_degree": int(QUADRATURE_DEGREE)}`), one compiled form per integrand — θ¹ − T∞ (both Δt), q_s, Πq_s, Π_L q_s — each value `allreduce`d before dividing by the allreduced ball volume. Print the ball volume against 4/3·π·a³.
  6. **Projection identity:** project θ¹ − T∞ (a V1 function, Δt₁, ccw) with `_projections` (tag `f"ops84_{arm}_id"`), and read |⟨Πθ⟩_B/⟨θ⟩_B − 1| at both tips.

**Gates (asserted; bands imported, none set).**
- **(A0m)** as I41: arm M global cells = `KEPT_ROD_CELL_RECORD`, tag-403 = `ROD_CELLS`, tag-3 equal to arm R's.
- **(B)** ≤ `B_BAND` on all four (arm, source) pairs.
- **(E)** `max_balance_residual` ≤ `E_BAND` on all eight solves.
- **Global (S2)** ≤ `S2_BAND` on all eight solves (predicted green: the projection keeps the domain mean, as entry 55 and U127 read).
- **Reproduction control:** the ccw Δt₁ global ratio ρc·θ̄¹/(Δt₁·G·P/V) equals `OPS-81` W2's printed 0.999994872614 (arm M, `20261105T183330Z_OPS-81-W2-i41.log:3683`) and 0.999995000025 (arm R, `:3728`) **within 1e-11** — ten units of the `.12f` print's last digit, the bound named here per the 18:00 ruling on I41. A miss means this route is not I41's: stop, BLOCKED.
- **Projection identity** ≤ 1e-10 at both tips on both arms (`OPS-83`'s LU round-off ceiling, `PROJECTION_ID_TOL`).
- max\|Im θ\| ≤ 1e-12 printed beside (I39 prints it; not a new gate).

**Printed (the separating measurement; no band set).** For each arm × source × tip, at Δt₁ and Δt₁/10: R_DG0 = \|ρc⟨θ¹⟩/(Δt⟨q⟩) − 1\| (the gate as written), R_Π against ⟨Πq⟩, R_L against ⟨Π_L q⟩, each Δt-ratio R(Δt₁)/R(Δt₁/10), and ⟨Πq⟩/⟨q⟩ − 1. **The pre-stated reading**, one line per arm × source × tip, printed, never asserted:
- **"form ACCOUNTS"** — R_DG0 > `S2_BAND` with Δt-ratio in `CONST_RATIO` [0.5, 2], and R_Π ≤ `S2_BAND`;
- **"form + O(Δt)"** — R_DG0 > `S2_BAND` Δt-independent, R_Π > `S2_BAND` with Δt-ratio ≥ `SCALING_RATIO` (5): the remainder is the step's conduction term (ruling 5), which a smaller Δt₁ would bring under the band;
- **"RESIDUAL"** — R_Π > `S2_BAND` with Δt-ratio < 5: something Δt-independent beyond the projection;
- **"INCONCLUSIVE"** — R_DG0 ≤ `S2_BAND` (this tip's source is too smooth to show the gap), or anything else.
Predicted (unlogged review estimate, not a band): R_DG0 between 1e-4 and 1e-2 on both arms — U127's and the gel twin's range — and the readings "form ACCOUNTS" or "form + O(Δt)".

**Window.** First the pre-check, a bare call: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. W0 (coding tests, knob unset): `scripts/testing/run_and_log.sh OPS-84-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_ops84_s2tip_form_rod_pair.py -v -s --tb=short'"`. W1: `scripts/testing/run_and_log.sh OPS-84-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_OPS84=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ops84_s2tip_form_rod_pair.py -v -s --tb=short'"` — `OPS-81` W2's form (`20261105T183330Z_OPS-81-W2-i41.log:12`) with this module. Both bare and foreground, Bash timeout 660000 ms. **Price:** I41 ran arms M 75.4 / R 66.0 / control 64.4 s, 306.0 s total (`…OPS-81-W2-i41.log:3772`); without the control arm and with six extra BE steps and five projections per arm on ≈ 52–67 k-cell submeshes, ≈ 250–330 s. If W1 is killed at 600 s, it is a size kill: park, never raise the timeout (drop the P1 source first on a re-pose).

**Tier.** Heavy, `-n 2`, complex. **Class:** lane-1 revalidation (diagnosis). **Cost:** ≈ 3.9 measured (the lane-1 median, 10:30 journal §6) / 25 predicted min. **Family:** `OPS-84`, attempt 1 of 4. It is not an attempt on U127's, entry 55's or `IMP-3`'s families: it executes no step and moves no gate.

**Traps already paid for.**
- `ArityMismatch` is a `BaseException` raised on rank 0 alone, so a raising rank hangs its peer to the ceiling (item 849). Conjugate every test function (`ufl.inner(u, v)`), and never put a deliberate refusal inside `pytest.raises` around `fem.form`.
- Comparisons on complex operands go through `ufl.real` (`OPS-22`); the ball indicator already does.
- A `SpatialCoordinate`-bearing integral without a pinned `quadrature_degree` can send FFCx into a compile that does not finish (`POST-5` step 1); pin it on every ball form.
- `_projections` builds a `LinearProblem` with a fixed options prefix per tag; reuse of a tag across calls collides.
- `_arm_with_403` monkeypatches `_rb.build_four_port_sweep` and restores it in `finally`; do not wrap it a second time. It also sets `TWO_TAG_ACK_KEY`, which `OPS-80`'s guard needs for the kept tag 403 (`OPS-81`).
- Reduce before asserting (`allreduce` on every ball integral and volume); count pytest verdicts per rank by token; use `-s`.
- Two heavy windows never run at once: `pgrep -c python3` reads 0 first. Harness calls are bare, foreground, `timeout -k 30`.

**Scope.** One mesh pair (I41's, ragged-sheet W1 record family, degree-1 EM, P1 thermal), two sources, one BE step. **Not claimed:** a re-pose of U127 or entry 55's gate (the weekly's), a band, a reading on U127's exact-rectangle mesh or on F-clinical, a fix to `_thermal_failures`, or any statement about I41's θ ratio.

**Negative result.** A red on (A0m), (B), (E), the global (S2), the reproduction control or the projection identity means the route is not I41's: print everything, stop, BLOCKED, journal. Any of the four readings is a result, not a failure: write it into the known-issues entry's Cause line — the "Not measured" clause becomes measured on a field-derived source, or the form reading is refuted there — in the same commit. No `PREDICTED RED:` line.

**Status it can move.** The known-issues entry of 2026-11-06: its Cause line's open clause is answered on a degree-1 \|E\|² beside a PEC rod tip and its gel twin. The 11-07 weekly's re-pose of the tip gate (U127, entry 55) gets field-derived evidence — whether a projected comparand alone, or a projected comparand at a smaller Δt₁, brings the tip under `S2_BAND`. `OPS-84` ✅.

**Depends on:** nothing.

**Outcome (12:00 CST slot, 2026-11-06).** ✅ `20261106T180958Z_OPS-84-W1.log`: 14 passed in 264.9 s at `-n 2` (harness 267 s; arms M 89.1 / R 75.1 s, total 263.2 s, `:3763`), heavy tier, complex build; W0 `20261106T180947Z_OPS-84-W0.log` (3 passed, 7 skipped, 3 s). (A0m) 182 371 / 15 000 / 51 693 = arm R's 51 693 (`:3645`); (B) ≤ 4.1e-14; (E) ≤ 3.8e-10; global (S2) ≤ 5.22e-06 on all eight solves; reproduction control — arm M 0.999994872611768 vs 0.999994872614 (2.2e-12), arm R 0.999995000024765 vs 0.999995000025 (2.3e-13) (`:3744–3745`); projection identity ≤ 4.4e-16 at both tips on both arms; max\|Im θ\| 0. Ball volumes 0.999929 / 1.000993 (M) and 0.981846 / 0.982903 (R, the rod hole) × 4/3·π·a³. Separating table (`:3747–3754`), Δt₁ → Δt₁/10 with Δt-ratio:

| arm src tip | R_DG0 | R_Π | R_L | ⟨Πq⟩/⟨q⟩ − 1 | reading |
|---|---|---|---|---|---|
| M ccw top | 9.4354e-03 → 9.5708e-03 (0.986) | 1.5196e-04 → 1.5265e-05 (9.954) | 1.0148e-01 (0.999) | −9.5859e-03 | form + O(Δt) |
| M ccw bottom | 1.6451e-04 → 7.3271e-05 (2.245) | 1.0142e-04 → 1.0187e-05 (9.956) | 9.9662e-02 (0.999) | +6.3084e-05 | INCONCLUSIVE |
| M P1 top | 1.8371e-03 → 1.7105e-03 (1.074) | 1.4052e-04 → 1.4125e-05 (9.948) | 4.8508e-02 (0.998) | +1.6963e-03 | form + O(Δt) |
| M P1 bottom | 5.9023e-03 → 5.8461e-03 (1.010) | 6.2204e-05 → 6.2569e-06 (9.942) | 4.4469e-02 (0.999) | +5.8398e-03 | form ACCOUNTS |
| R ccw top | 4.0523e-03 → 3.9690e-03 (1.021) | 9.2880e-05 → 9.2323e-06 (10.060) | 2.7292e-02 (0.997) | −3.9598e-03 | form ACCOUNTS |
| R ccw bottom | 2.9397e-03 → 3.0233e-03 (0.972) | 9.2589e-05 → 9.2618e-06 (9.997) | 2.9794e-02 (0.997) | +3.0326e-03 | form ACCOUNTS |
| R P1 top | 1.3758e-03 → 1.2818e-03 (1.073) | 1.0453e-04 → 1.0398e-05 (10.053) | 4.7682e-02 (0.998) | −1.2714e-03 | form + O(Δt) |
| R P1 bottom | 3.7808e-03 → 3.8794e-03 (0.975) | 1.0914e-04 → 1.0918e-05 (9.996) | 5.0440e-02 (0.998) | +3.8904e-03 | form + O(Δt) |

Reading: R_DG0 sits in the predicted 1e-4 … 1e-2 range and equals the projection's ball-mean shift; once the comparand is ⟨Πq⟩ the remainder is O(Δt) (ratio ≈ 10 per decade) and under `S2_BAND` on every row at Δt₁/10. No row reads RESIDUAL. The one INCONCLUSIVE row (M ccw bottom) has R_DG0 just over the band with a non-constant Δt-ratio, because its form gap (6.3e-05) is smaller than the O(Δt) term. The lumped projection is worse than DG0, so it is not the comparand. Implementation choice, coded in the module's `_reading`: the reading uses R_Π at Δt₁ (the gate's step). Scope as stated: ragged-sheet I41 pair only; no band, no re-pose.

**Audit PASS — 2026-11-06 18:00 review (`auditor`).** Footer Status 0, 267 s; `Commit:` header `71c41460` = parent of `c18386ce` (`20261106T180958Z_OPS-84-W1.log:3`, `:3777–3780`). 14 passed per rank by token (`:3772`, `:3775`; 28 `PASSED` tokens, no FAILED / ERROR). Asserts executed at `test_ops84_s2tip_form_rod_pair.py:303–306` (A0m), `:314–317` ((B), (E), global (S2)), `:324` (reproduction), `:331` (projection identity); `ASSERTED` lines at `:3679–3688`, `:3725–3734`, `:3744–3745`. Every quoted digit traced (reproduction 2.232e-12 / 2.346e-13, `:3744–3745`; table `:3746–3754`; readings 3 / 4 / 1 / 0, `:3755–3762`). Bands imported unchanged; the one new constant, `REPRO_TOL` = 1e-11, is the brief's. No `src/` change. **Caveat corrected in the §7 row:** "R_DG0 Δt-independent" held on the seven classified rows only; M ccw bottom (the INCONCLUSIVE row) reads a Δt-ratio of 2.245 (`:3749`). Noted, not edited: `_reading` also requires the R_DG0 ratio in `CONST_RATIO` for "form + O(Δt)", slightly stricter than the brief's wording.
