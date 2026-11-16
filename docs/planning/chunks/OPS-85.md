# OPS-85 — an exercising node for the _tuned_quadrature_temperature gradient refusal

*History file for the §7 `OPS-85` row (`OPS-46`). Opened by the 2026-11-08 03:00 daily review, lane 1, from U125's disclosure (3) as the 11-07 weekly ratified it (`PROJECT_PLAN.md` §10, "Public surface ratified": "the `src` gradient refusal without an exercising node — the last is a lane-1 test to write, `OPS-85` below"). U125 is `WF-1` step 28, §9 item 830. Its Outcome (`docs/planning/chunks/WF-1.md:3590`) says: "(3) The `_tuned_quadrature_temperature` gradient refusal is in `src` but has no node that exercises it, because reaching it needs a tuned Study."*

**The refusal.** `src/fem_em_solver/study.py:3396–3397`, inside `_tuned_quadrature_temperature` (`:3378`):

```
if transient.initial_gradient_k_per_m is not None:  # `WF-1` step 28 (§9 item 830)
    raise ValueError("initial_gradient_k_per_m is defined on a ThermalBoxGeometry only")
```

The same message also stands at `:1659` (the EM sphere, `_study_temperature`) and `:2139` (`_study_transient_temperature`). The sphere's site is exercised by `tests/package/test_wf15_study_sphere_transient.py:127–129`. **No node matches the message at all**: a grep of `tests/` for "ThermalBoxGeometry only" returns nothing, and the sphere node matches only `"initial_gradient_k_per_m"`. A node on this site must therefore show that it is *this* site that raised.

**The route, read against the code (no run).** A tuned two-port Study needs no mesh to reach the check:
- `Study.temperature` is `_study_temperature` (`:1734`).
- Its first branch takes `excitation._temperature_port_study` (`:1625–1629`).
- `TwoPortQuadratureDrive._temperature_port_study` (`:3030–3035`) calls `_tuned_quadrature_temperature`.
- The function checks `cooling` (`:3392`), `transient` (`:3394`) and the gradient (`:3396`) **before** `_tuned_result(study, comm)` (`:3412`). That call is the first that solves (`:3242–3246`).

`tests/package/test_wf16_study_tuned_f_human.py`'s `_study()` (`:129–137`) builds the tuned copper F-human Study with no mesh, and its (K0) node already calls refusals on it before any build (`:174–182`).

### §9 item 879 — 2026-11-08 (lane 1, OPS-85 — U125's disclosure (3); smoke)

**Executes** the 11-07 weekly's lane-1 target "`OPS-85` — the `_tuned_quadrature_temperature` gradient refusal gets an exercising node (U125's disclosure (3)), smoke". Family OPS-85, attempt 1 of 4.

**Code.** New `tests/package/test_ops85_tuned_gradient_refusal.py` only. No `src/` edit, no landed module edited, and no envelope line: an envelope line would move the generated pages (Wgen / Wval / Wapi) and lift this item out of smoke. The weekly may add one.
- **Imports:**
  - `_study` from `tests.package.test_wf16_study_tuned_f_human`;
  - `RobinCooling` and `fem_em_solver.study` as `study_mod`, with `study_mod.TransientSpec`;
  - `RHO_C` from `tests.validation.test_mat12_transient_sar_heating_sphere`;
  - `H_GEL`, `K_GEL`, `T_INF` from `tests.validation.test_mat15_coil_driven_heating_phantom`, as `test_wf1_yaml_f_human_thermal.py:55–64` imports them.
- **The call:** `st.temperature(RobinCooling(K_GEL, H_GEL, T_INF), transient=TransientSpec(dt_s=1.0, n_steps=1, rho_c_j_per_m3_k=RHO_C, initial_gradient_k_per_m=g), balls={"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]})` on `st = _study()`.
- **Fakes.** Monkeypatch `study_mod._tuned_result` with a recorder that counts calls and raises a local `_Reached` exception. Wrap `MeshGenerator.birdcage_port_domain` in a counter that raises `_Tripwire`, as `test_geo38_f_clinical_3t_matrix.py:150–165` does. Patch the `fem_em_solver.study` module object, because `_tuned_quadrature_temperature` looks `_tuned_result` up as a module global at call time.

**Gates (asserted).**
- **(Ref) the anchor.** With g = 10.0, `pytest.raises(ValueError)` holds, with `str(exc) == "initial_gradient_k_per_m is defined on a ThermalBoxGeometry only"` exactly. **`excinfo.traceback[-1].name == "_tuned_quadrature_temperature"`**: the raising frame is this site, not `:1659` or `:2139`. The `_tuned_result` recorder count is 0, and the builder count is 0.
- **(Ctl) the accepted control.** With `initial_gradient_k_per_m=None`, the same call raises `_Reached`, never `ValueError`, and the recorder count is 1. So every argument check before the mesh accepts this spec, and the refusal in (Ref) is the gradient's alone. The builder count stays 0.
- **(Ctl-g0)** g = 0.0 is refused like 10.0. The check is `is not None`, not truthiness (`:3396`), so a zero gradient is never silently accepted.
- **Both ranks** reach each verdict. The raises come before any collective, so there is no `MPI_Bcast` to hang on. Count the verdicts per rank, by token.

**Printed.** The message, the raising frame's name and line number (`excinfo.traceback[-1].lineno + 1`; pytest's is 0-based, printed only), and both counters.

**Windows.** Bare, foreground, Bash timeout 660000 ms. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 first.
- **W1:** `scripts/testing/run_and_log.sh OPS-85-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_ops85_tuned_gradient_refusal.py tests/package/test_wf16_study_tuned_f_human.py::test_wf16_s9a_api_and_refusals -v -s --tb=short > /workspace/logs/ops85-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ops85-w1-raw.log; cat /workspace/logs/ops85-w1-raw.log; exit \$rc'"`. This uses the durable capture of `20261104T030134Z_IMP-8-I55-W2-item803.log:12`, escaped as `docs/planning/chunks/GEO-38.md:233` spells it. The shape is `WF-15` step 16's no-solve refusal window (`20261107T180132Z_WF-15-s16-W0.log:12`, complex, `-n 2`, 5 passed / 2 skipped per rank, 5 s, `:70`, `:77`). The `_study()` source's own (K0) node rides along as the control that the imported helper still builds without a mesh.
- **W-real** (optional; run it if W1 is green before minute 30): the same command without `source /usr/local/bin/dolfinx-complex-mode` / `FEM_EM_REQUIRE_COMPLEX=1`, labelled `OPS-85-W-real`, raw log `ops85-wreal-raw.log`. The refusal is a plain-Python check in a build-agnostic function, so it should be green in both builds.

**Tier.** Smoke, `-n 2`, complex (W-real real), no mesh, no solve. **Class:** lane 1. **Cost:** ≈ 3.9 measured / 15 predicted min (W1 ≈ 5–10 s).

**Traps already paid for.**
- A complex-build refusal is never asserted through `fem.form` inside `pytest.raises`: `ArityMismatch` is a `BaseException`, raised on rank 0 alone, and the peer hangs to the ceiling (item 849, §9 rulings). This refusal is raised before any form, so it is safe, but nothing in the node may build a form.
- `-k a or b` splits into stray argv inside the quoted container command. Select by file path and `::node`, as W1 does.
- Pytest captures prints without `-s`.
- `.eval` is not involved.
- Do not monkeypatch through a re-imported copy of `study.py`. Patch `fem_em_solver.study`, the object `_tuned_quadrature_temperature`'s globals are.
- A headless session must not background a harness run.

**Scope.** One refusal site exercised, with its accepted control. **Not claimed:**
- a gradient on a tuned Study, or any tuned thermal solve;
- the sites at `:1659` / `:2139`;
- an envelope or API-page line;
- U125's other two disclosures, which are refusals already pinned by item 830's W0.

**Negative result.**
- **(Ref) red with the message raised from another frame** means the tuned route reaches a different check first: a finding about the route. Print the frame, stop, BLOCKED.
- **(Ref) red because nothing raises** means the refusal is unreachable as the disclosure assumed. File it in known-issues; never edit `study.py` here.
- **(Ctl) red, with `ValueError` from an earlier check,** means the control spec is wrong for the route. Fix the test's arguments, keep the first log, and re-run W1 once (rule (i)).
- No `PREDICTED RED:` line.

**Status it can move.**
- `OPS-85` ⬜ → ✅: a refusal pinned by an exact-message, exact-frame assertion with an accepted control and zero builder calls. The closure's audit rules whether that meets §4's "quantitative assertion" clause, since the gates are exact identities, not a closed form. If it does not, the row is 🟡 with the node landed.
- U125's disclosure (3) is retired in the `WF-1` row's next touch.

**Depends on:** nothing.

**Outcome (09:00 slot, 2026-11-08): complete, ✅ pending the audit's §4 ruling.** New module `tests/package/test_ops85_tuned_gradient_refusal.py` only (no `src/`, envelope or landed-module edit). (Ref) g = 10.0 and (Ctl-g0) g = 0.0 each raise `ValueError` with the exact message from the frame `_tuned_quadrature_temperature` (printed line 3534 — the site moved from `:3396` with item 878's `study.py` edit), `_tuned_result` calls 0, builder calls 0, on both ranks. (Ctl) g = None raises `_Reached` from the recorder, calls 1, builder 0. W1 (complex, `-n 2`): `20261108T152052Z_OPS-85-W1.log`, 8 passed per rank (`test_complex_mode.py` 4, this module 3, WF-16's (K0) 1), 3 s, `[capture] rc=0`. W-real (real, `-n 2`): `20261108T152105Z_OPS-85-W-real.log`, 5 passed / 3 skipped per rank (the three complex-only environment nodes), 3 s. U125's disclosure (3) can be retired at the `WF-1` row's next touch.

**Audit (2026-11-08 10:30 daily review): PASS.** The `auditor` traced every claim to `20261108T152052Z_OPS-85-W1.log` (per-rank verdict tokens `:63`, `:65`, `:66`, `:69`, `:70`, `:72`; 8 passed twice `:79`, `:82`; `[capture] rc=0` `:83`; Status 0, 3 s `:85–88`) and `20261108T152105Z_OPS-85-W-real.log` (5 passed / 3 skipped twice `:77`, `:80`; Status 0 `:83–86`). The asserts are `test_ops85_tuned_gradient_refusal.py:77–96`, and the refusal is `study.py:3533–3534`. **§4 ruling:** exact string equality on the message, exact equality of the raising frame's name and integer call-count identities (0 / 0 on the refusals, 1 / 0 on the control) are not the finiteness-only kind §4 forbids. Each is falsifiable: delete the refusal, change the message, move it, or make it truthy, and an assert fails. They meet §4's quantitative clause in its operative sense, though not by a literal match to its bullet list (closed form, rate, conservation / reciprocity identity). No earlier closure rested on a refusal alone (`OPS-81` paired its refusal contract with reproduced records), so **this sets a precedent, which the 11-09 weekly ratifies or reverts**; if reverted, the row returns to 🟡. Caveat: builder calls 0 follows from `_tuned_result` calls 0 under the stub and adds no independent evidence.
