# OPS-83 — the (S2-tip) comparand's form on a rough DG0 source: a stub measurement that separates the two readings

*History file for the §7 `OPS-83` row (`OPS-46`). Opened by the 2026-11-06 03:00 daily review (lane 1, known-issues diagnosis) from the known-issues entry of the same date: `xl` entry 55's (S2-tip) red, which `WF-15` step 14 / U127 met first. The weekly re-poses the gate; this chunk only measures.*

**Why.** Two windows have failed the same assertion on field-derived sources, and the gel twin of entry 55 misses it too. The global (S2) on those runs is green. The gate is |ρc·⟨θ¹⟩_B/(Δt₁·⟨q⟩_B) − 1| ≤ `S2_BAND`. One reading is that the comparand has the wrong form: one BE step from rest gives ρc·θ¹/Δt₁ → M⁻¹b as Δt → 0, the consistent-mass L2 projection Πq of the DG0 q onto P1, and Πq keeps the domain mean but not a ball mean. The other reading is physics the gate is meant to catch. Under the form reading, (i) the gap needs a rough source, (ii) it does not scale with Δt₁, and (iii) it vanishes against ⟨Πq⟩_B. U127's Δt₁/10 control already saw (ii) (`log-pathologist`, 2026-11-06). The pathologist named this measurement as the cheap check.

**Done-when.** `tests/validation/test_ops83_s2tip_form_on_stub.py` lands, with every asserted gate below green at `-n 2` in the complex build, and the separating table printed. The known-issues entry's Cause line is updated with the measured reading, confirmed or refuted, and the entry stays open; it resolves with the weekly's re-pose ⇒ ✅.

### §9 item 859 — 2026-11-06 (lane 1, known-issues diagnosis; implementer; standard, `-n 2`)

**Code.** New `tests/validation/test_ops83_s2tip_form_on_stub.py`, complex build. No `src/` edit; no landed module edited.
- **Imports, unchanged, from `test_geo36_rod_f_clinical_heating`:** `_columns_e_sup`, `_merge_rod_into_phantom`, `_tag_submesh`, `_outer_robin_tags`, `_heat_source`, `_transient`, `_ball_reader`, `DT_1_BOX`, `STUB_CENTRE`, `STUB_PROBES`, `RHO_C`, `T_INF`, `S2_BAND`, `E_BAND`, `QUADRATURE_DEGREE`, `ROBIN_TAG`, `PHANTOM_CELL_TAG`, `ONE_GRAM_KG`, `PHANTOM_RHO_KG_PER_M3`, `averaging_ball_radius`, `_allsum`, `_say`. Also `_k9_fixture` from `test_geo36_rod_f_clinical_drive`.
- **Read every signature first.** Copy K10es's set-up exactly (`test_geo36_rod_f_clinical_heating.py:656–666`): the stub mesh, `_columns_e_sup` under `monkeypatch.context()`, the merged arm, `_tag_submesh` and `_outer_robin_tags` with `STUB_CENTRE`. Do not call `_thermal_arm`: its (B) compares against G·P, which a hand-made source does not satisfy.
- **Sources** (DG0 on the submesh):
  - **(u)** q_u = `_heat_source(…)` with G = 1/P_gel, as K10es does.
  - **(r)** q_r = q_u·(1 + ½·s), where s = sign(sin(πx/h₀)·sin(πy/h₀)·sin(πz/h₀)) at each cell midpoint and h₀ = 0.1/12 m, the stub's `create_box` cube size (`test_geo36_rod_f_clinical_drive.py:1188`). Every tet of one hex cube gets the same sign, which gives a rank-independent checkerboard with cell-scale jumps of the kind a degree-1 \|E\|² has. Print its owned min / max, allreduced.
- **Projection:** Πq solves (u, v) = (q, v) on the temperature space V1 = `res.temperature.function_space`. Use `ufl.inner(u, v)` and `ufl.inner(q, v)` (conjugated test function, complex build) and a direct LU (`petsc_options_prefix` unique). Also form the lumped projection Π_L q = b / (M·1), printed only.
- **Runs:** BE one step from T∞ with `_transient(sub, ft, q, dt=DT_1_BOX, n_steps=1, theta=1.0)`, for q ∈ {q_u, q_r} and Δt ∈ {`DT_1_BOX`, `DT_1_BOX`/10}. That is four solves on a 12³ stub.
- **Ball reads:** `_ball_reader(sub, {"th1": θ¹ − T∞, "q": q, "pq": Πq, "plq": Π_L q})` at the stub tip `fx[4]["a"]` (K10es's `tips`) and at `STUB_PROBES["interior"]`, radius a_1g. If a ball is not wholly in the submesh, print its volume against 4/3·π·a³.

**Gates (asserted; bands imported).**
- **(E)** `max_balance_residual` ≤ `E_BAND` on all four solves.
- **(S2) global**: |ρc·θ̄¹/(Δt·q̄) − 1| ≤ `S2_BAND` on all four solves, with θ̄¹ and q̄ as submesh means. The projection keeps the domain mean, so this is predicted green on (r), as it read green in entry 55 and U127.
- **Reproduction control:** on (u) at `DT_1_BOX`, the tip value of |ρc⟨θ¹⟩_B/(Δt⟨q⟩_B) − 1| equals K10-S's printed tip 8.0400e-06 (`20261101T103929Z_GEO-36-I32a-W0b.log:76`) within 5e-11, its last printed digit.
- **Projection identity:** on (u), |⟨Πq_u⟩_B/⟨q_u⟩_B − 1| ≤ 1e-10 at both balls. A field in V1 projects to itself, and q_u is uniform on the stub (K10-q). The bound is the LU round-off ceiling, not a physics band.

**Printed (the separating measurement; no band set).** For each source, Δt and ball:
- R_DG0 = |ρc⟨θ¹⟩/(Δt⟨q⟩) − 1|, which is the gate as written;
- R_Π = |ρc⟨θ¹⟩/(Δt⟨Πq⟩) − 1|;
- R_L = |ρc⟨θ¹⟩/(Δt⟨Π_L q⟩) − 1|;
- the Δt-ratio R(Δt₁)/R(Δt₁/10) for each.
**The pre-stated reading** (printed as one line, never asserted):
- **"form CONFIRMED"** if on (r) R_DG0 > `S2_BAND` with Δt-ratio in [0.5, 2] (Δt-independent), and R_Π ≤ `S2_BAND` with Δt-ratio ≥ 5 (an O(Δt) residual);
- **"form REFUTED"** if R_Π > `S2_BAND` on (r);
- **"INCONCLUSIVE"** otherwise, including when R_DG0 on (r) is itself ≤ `S2_BAND` (the stub is too smooth to show the gap).
The predicted values come from unlogged review arithmetic. R_DG0(r) should be at the per-cent level, because the ball (a_1g ≈ 6.2 mm) is smaller than the 8.3 mm cube and straddles checkerboard jumps. R_Π(r) should be at the 1e-5 level of K10-S's slope (2.0247e-05).

**Window.** W1: `scripts/testing/run_and_log.sh OPS-83-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ops83_s2tip_form_on_stub.py -v -s --tb=short'"`. This is a bare foreground call with Bash timeout 660000 ms. The form is `20261101T103929Z_GEO-36-I32a-W0b.log:12` with this module added. **Price:** I32a's whole stub module ran 7 tests in 14 s at `-n 2` (`…W0b.log:133`). Four BE solves plus two projections should take ≈ 20–60 s.

**Tier.** Standard, `-n 2`, complex. **Class:** lane-1 revalidation (diagnosis). **Cost:** ≈ 6.0 measured / 20 predicted min. **Family:** `OPS-83`, attempt 1 of 4. It is not an attempt on U127's or entry 55's family: it executes neither step and moves neither gate.

**Traps already paid for.**
- `ArityMismatch` is a `BaseException` and `mpi_jit` raises it on rank 0 alone, so a raising rank hangs its peer to the ceiling (item 849). Conjugate every test function and never put a deliberate complex refusal inside `pytest.raises` around `fem.form`.
- `ufl.lt` on complex operands is wrapped in `ufl.real` in `_ball_reader` already; any new comparison must do the same (`OPS-22`).
- `sign(sin…)` is computed in numpy on the midpoints, never in UFL. Midpoints come from `dolfinx.mesh.compute_midpoints` over owned + ghost cells, and the DG0 dofs are set from them, then `scatter_forward`.
- Reduce before asserting; count verdicts per rank by token; use `-s`.
- K10es runs under `monkeypatch.context()` because `_columns_e_sup` rebinds module names. Do the same, and do the ball reads inside the context.

**Scope.** A stub measurement on a 12³ box. **Not claimed:** a re-pose of U127 or of entry 55's gate (the weekly's), a new band, any clinical-scale number, or a fix to `_thermal_failures`.

**Negative result.** A red on (E), the global (S2), the reproduction control or the projection identity means the stub route is not K10es's: print everything, stop, BLOCKED. Any of the three readings is a result, not a failure: write it into the known-issues entry's Cause line (hypothesis → measured, or refuted) in the same commit. No `PREDICTED RED:` line.

**Status it can move.** The known-issues entry of 2026-11-06 gains a measured Cause, or a refuted hypothesis. The 11-07 weekly's re-pose of the tip gate (U127, entry 55, and entry 61's EM twin if it ever gains heating) gets its evidence. `OPS-83` ✅.

**Depends on:** nothing.

**Outcome (2026-11-06, 06:00 CST slot — 12:00 UTC; label corrected by the 10:30 review): DONE, form CONFIRMED.** *Audited PASS by the 2026-11-06 10:30 review (`20261106T120232Z_OPS-83-W1.log:68–84`, `:92–95`).* W1 `20261106T120232Z_OPS-83-W1.log` (standard, `-n 2`, complex; Status 0, 13 s; 5 passed per rank). Gates: (E) ≤ 3.450e-10; global (S2) 2.0247e-05 (Δt₁) / 2.0249e-06 and 2.0246e-06 (Δt₁/10) on both sources; K10-S tip 8.039988e-06 vs 8.0400e-06 (|diff| 1.183e-11); projection identity 6.661e-16 at both balls. Table (tip / interior; Δt₁ → Δt₁/10, ratio): (u) all three R coincide, 8.0400e-06 → 8.0980e-07 (9.928) / 8.3513e-06 → 8.3959e-07 (9.947). (r) R_DG0 1.5793e-03 → 1.5732e-03 (1.004) / 2.1076e-03 → 2.1146e-03 (0.997); R_Π 6.8172e-06 → 6.8552e-07 (9.944) / 7.9282e-06 → 8.5219e-07 (9.303); R_L 5.1545e-05 → 4.5413e-05 (1.135) / 1.5441e-05 → 8.3648e-06 (1.846). q_r on [500, 1500] W/m³; ball volumes 1.005418 / 0.999642 × 4/3·π·a³ (both balls inside the submesh). The predicted R_DG0(r) "per-cent level" read at 0.16–0.21 %; R_Π(r) at the 1e-5 level as predicted. Not claimed: that the clinical-scale reds are this gap alone; a re-pose; a band.
