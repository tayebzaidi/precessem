"""`WF-17` step 43 (§9 item 941, §10 U152): the implant-rod page's one-step study file.

``docs/guide/implant_rod_in_a_study.md`` "One step from the study file" carries an
``implant-rod-heating-step`` block (``run_study_file`` on the tracked
``tests/package/wf15/rod_heating.yaml`` with ``source_scale = 1/P`` and its
``transient:`` replaced by one backward-Euler step of Δt₁/10) and an
``implant-rod-heating-check`` block (``femem check`` on the tracked
``tests/package/wf15/malformed_rod_heating_duty.yaml``).

Brief: docs/planning/chunks/WF-17.md, "§9 item 965" (§10 U152′, re-posed by the
2026-11-11 weekly's ruling, option (d), which binds over item 941's brief; and the
2026-11-10 10:30 / 18:00 review rulings on records vs readings).

Gates (bands imported, none set):
- W0, no mesh: (X) one step block, one ``run_study_file(``, no ``Study(`` / ``probe``,
  no ``yaml`` fence beyond ``implant-rod-file``; (In) the prose's four transient
  values equal ``_be(DT_1 / 10.0)``'s fields; (Ref) the malformed file is refused at
  ``thermal.duty``, ``femem check`` rc 2 on both ranks; (Ref-Py) ``Study.temperature``
  with ``duty=`` on the rod Study raises ``ValueError`` before any mesh; (D) 0 dead
  docrefs; (c-page) a removed or doubled step / check tag raises.
- W1 (``FEM_EM_WF17_S43=1``, complex, ``-n 2``): (Scale) / (Scale CN) preconditions;
  (Cells) ×2; (E) ×4; (Y-θ) per tip, the CN study file's θ(t*) at the tip balls equal
  to the ``Study``'s at ``IDENTITY_BAND`` (item 908's gate, re-run here); (Band) R_Π
  per tip and |s − 1| ≤ ``S2_BAND`` on both routes, with r = ρc⟨θ¹⟩_B/(Δt⟨Πq⟩_B)
  against one comparand; ``R_Pi form``.  Controls (asserted as misses): (c-dt) the
  YAML step's r against the Study's Δt₁ step; (c-θ) the CN file's top-tip θ(t*)
  against the Study's bottom-tip θ(t*).
- Printed with no band (the 11-10 18:00 ruling: a 1e-12 slope identity is not
  re-posed): (Y-Π) |r_Y − r_S| per tip and (Y-S2) |s_Y − s_S|, each beside the same
  route pair's (E); the volume-mean θ(t*) on the Study and the CN file.
"""

from __future__ import annotations

import os
import re
import time

import pytest
import yaml
from mpi4py import MPI

import tests.package.test_wf17_implant_rod_tutorial as _rod_tut
from tests.complex_mode import IS_COMPLEX
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf15_study_implant_rod import ROD_EXACT_CELL_RECORD
from tests.package.test_wf15_study_rod_heating import _ball_comparands, _be, _yaml_mapping
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf17_implant_rod_heating_tutorial import HEAT_FILE_TAG, HEAT_TAG, _tips
from tests.package.test_wf17_implant_rod_tutorial import (
    MIN_GUIDE_PAGES,
    PAGE,
    _exec_block,
    _scratch,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_geo36_rod_f_clinical_heating import _tip_slope_rows
from tests.validation.test_imp3_rod_birdcage import TIPS
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    DT_1,
    E_BAND,
    RHO_C,
    S2_BAND,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

S43_ENV = "FEM_EM_WF17_S43"
STEP_TAG = "implant-rod-heating-step"
CHECK_TAG = "implant-rod-heating-check"
STEP_FILE = "rod_heating_step.yaml"
SECTION = "## One step from the study file"
REPO = PAGE.resolve().parents[2]
MALFORMED = REPO / "tests" / "package" / "wf15" / "malformed_rod_heating_duty.yaml"
SCRATCH = REPO / "logs" / "wf17" / "s43"
DRIVE = "P1"
#: Item 884's printed priors (`20261108T223537Z_WF-15-s14p-W1.log:5674–5678`): printed
#: beside the readings, never asserted here.
PRIOR_P = "1.646387247816e-06"
PRIOR_RPI_10 = {"top": "1.3394e-05", "bottom": "9.6264e-06"}
PRIOR_S2_10 = "4.8365e-07"
#: Item 908's (Y) rel on the same pair (`20261109T103942Z_WF-17-s32p-W1.log:3725`, `:3727`),
#: printed beside the (Y-θ) readings, never asserted here.
PRIOR_Y_THETA = {"top": "2.887e-15", "bottom": "6.661e-16"}
ROD_YAML = REPO / "tests" / "package" / "wf15" / "rod_heating.yaml"
DUTY_MSG = "duty / window_s / sense are not defined on the rod four-port Study"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S43_ENV, "").strip()
    if raw not in ("", "1"):
        raise ValueError(f"{S43_ENV}={raw!r}: only unset/empty or '1' are accepted")
    return raw == "1"


def _section(text: str) -> str:
    parts = text.split(SECTION)
    if len(parts) != 2:
        raise LookupError(f"{len(parts) - 1} sections {SECTION!r}; need exactly one")
    return parts[1].split("\n## ")[0]


def _prose(section: str) -> str:
    """The section with fenced blocks and link targets removed."""
    no_fence = re.sub(r"```.*?```", "", section, flags=re.S)
    return re.sub(r"\]\([^)]*\)", "]", no_fence)


def _prose_transient(section: str) -> dict:
    vals = dict(re.findall(r"`(dt_s|n_steps|theta|startup_half_steps): ([^`]+)`", _prose(section)))
    assert set(vals) == {"dt_s", "n_steps", "theta", "startup_half_steps"}, vals
    return vals


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s43_step_block_and_refusal(capsys):
    from fem_em_solver import cli
    from fem_em_solver.study_yaml import StudyConfigError, load_study

    comm = MPI.COMM_WORLD
    text = PAGE.read_text()
    block = extract_tagged_block(text, STEP_TAG)
    check = extract_tagged_block(text, CHECK_TAG, lang="bash")
    section = _section(text)

    # (X) one step block: one run_study_file(, no Study(, no probe; the file name once.
    assert block.count("run_study_file(") == 1, block
    assert "Study(" not in block and "probe" not in block, block
    assert block.count(f'"{STEP_FILE}"') == 1 and text.count(f'"{STEP_FILE}"') == 1
    assert check.strip() == f"femem check {MALFORMED.name}", check
    yaml_fences = [ln for ln in text.splitlines() if ln.startswith("```yaml")]
    assert yaml_fences == ["```yaml implant-rod-file"], yaml_fences
    assert "```yaml" not in section
    assert f"](../../tests/package/wf15/{MALFORMED.name})" in section
    assert "thermal.duty" in _prose(section)

    # (In) the prose's four transient values = _be(DT_1 / 10.0); no other digit in the prose.
    # Measured (W0 r1, `20261110T210357Z_WF-17-s43-W0.log:72`): the brief's `dt_s: 3.135`
    # read ('3.135', 0.0015675), rel 1999. DT_1 is the calorimetric step 2.5e-6 T_R
    # = 0.015675 s, not the CN step 31.35 s, so the page states 0.0015675; the gate is unchanged.
    be10 = _be(DT_1 / 10.0)
    vals = _prose_transient(section)
    dt_rel = abs(float(vals["dt_s"]) / be10.dt_s - 1.0)
    assert dt_rel <= IDENTITY_BAND, (vals["dt_s"], be10.dt_s)
    assert int(vals["n_steps"]) == be10.n_steps
    assert float(vals["theta"]) == be10.theta
    assert int(vals["startup_half_steps"]) == be10.startup_half_steps
    digits = re.findall(r"\d+(?:\.\d+)?", _prose(section))
    assert digits == [vals["dt_s"], vals["n_steps"], vals["theta"], vals["startup_half_steps"]], digits
    with capsys.disabled():
        _p0(f"\n[WF-17 s43] (X) step block one run_study_file, no Study(/probe; yaml fences "
            f"{yaml_fences}")
        _p0(f"[WF-17 s43] (In) prose {vals} vs _be(DT_1/10) dt {be10.dt_s!r} (rel {dt_rel:.3e}), "
            f"n {be10.n_steps}, theta {be10.theta!r}, half {be10.startup_half_steps}; prose "
            f"digits {digits}")

    # (Ref) the malformed file: thermal.duty refused on load and by `femem check` (rc 2).
    with pytest.raises(StudyConfigError) as err:
        load_study(MALFORMED)
    field_ok = err.value.field == "thermal.duty"
    capsys.readouterr()
    rc = cli.main(["check", str(MALFORMED)])
    check_out = capsys.readouterr()
    rc_all = (comm.allreduce(int(rc), op=MPI.MIN), comm.allreduce(int(rc), op=MPI.MAX))
    with capsys.disabled():
        _p0(f"[WF-17 s43] (Ref) load_study field {err.value.field!r}, reason {err.value.reason!r}; "
            f"femem check rc (min, max) {rc_all}; check stdout {check_out.out!r} stderr "
            f"{check_out.err!r}")
    assert field_ok, err.value.field
    assert rc_all == (2, 2), rc_all

    # (N) the page's limits stay; no external tool named.
    for token in ("No two-port drive", "No 1 g rung", "No published comparand", "No tissue"):
        assert token in text, token
    for token in ("AED", "HFSS"):
        assert token not in text, token

    # (c-page) each new tag removed or doubled raises.
    for tag, lang, body in ((STEP_TAG, "python", block), (CHECK_TAG, "bash", check)):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag, lang=lang)
        with pytest.raises(LookupError):
            extract_tagged_block(text + f"\n```{lang} {tag}\n{body}```\n", tag, lang=lang)
    _p0("[WF-17 s43] (c-page) removed / doubled step and check tags raise LookupError")


def test_wf17_s43_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s43] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


def test_wf17_s43_rod_route_refuses_duty_in_python():
    """(Ref-Py) the rod Study refuses ``duty=`` in Python before any mesh (`study.py:1650–1651`).

    Plain Python logic, rank-uniform (no ``fem.form``), so ``pytest.raises`` is safe.
    """
    from fem_em_solver import DutyCycle
    from fem_em_solver.study_yaml import _load

    study, _kind, _outputs, rod = _load(ROD_YAML)
    with pytest.raises(ValueError, match=re.escape(DUTY_MSG)) as err:
        study.temperature(rod.cooling, transient=rod.transient, balls=rod.balls, drive=rod.drive,
                          duty=DutyCycle(period_s=250.8, duty=0.5))
    unset = study.result is None and getattr(study, "fields_result", None) is None
    _p0(f"\n[WF-17 s43] (Ref-Py) {type(err.value).__name__}: {str(err.value)!r}; result / "
        f"fields_result unset afterwards: {unset} (ASSERTED)")
    assert unset, (study.result, getattr(study, "fields_result", None))


# ------------------------------------------------------------- solve window (-n 2)


def _ratios(theta1: dict, dt: float, comps: dict) -> dict:
    """r = ρc⟨θ¹⟩_B / (Δt⟨Πq⟩_B) per tip; ``_tip_slope_rows``'s R_Π is |r − 1|."""
    return {n: RHO_C * th / (dt * comps[n]["pq"]) for n, th in theta1.items()}


def test_wf17_s43_step_block(monkeypatch):
    if not _knob():
        pytest.skip(f"{S43_ENV} unset: the one-step file's solve window is environment-selected")
    import numpy as np
    import ufl
    from dolfinx import default_scalar_type, fem

    import fem_em_solver.study as study_module

    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port Study needs the complex build"
    assert comm.size == 2, comm.size
    monkeypatch.setattr(_rod_tut, "SCRATCH_ROOT", SCRATCH)
    text = PAGE.read_text()
    names = list(TIPS)

    # 1. The page's heating block: one EM solve, P, cooling (scale 1/P), study.
    t0 = time.perf_counter()
    ns, _out = _exec_block(monkeypatch, extract_tagged_block(text, HEAT_TAG), "heating")
    t_py = time.perf_counter() - t0
    st, cooling, heat = ns["study"], ns["cooling"], ns["heat"]

    # 2. Study side: one BE step at DT_1/10 (the gate) and at DT_1 (the control).
    t0 = time.perf_counter()
    balls = {"mass_g": 1.0, "centres": ns["tips"]}
    be10 = st.temperature(cooling, transient=_be(DT_1 / 10.0), balls=balls, drive=DRIVE, comm=comm)
    be = st.temperature(cooling, transient=_be(), balls=balls, drive=DRIVE, comm=comm)
    t_be = time.perf_counter() - t0

    # 3. One comparand for both routes; V_sub from the Study's phantom submesh.
    t0 = time.perf_counter()
    res = st.fields_result
    comps = _ball_comparands({"res": res, "G": 1.0 / ns["P"], "be10": be10}, None)
    sub, _pid = study_module._tag_submesh(res.mesh, res.cell_tags, PHANTOM_CELL_TAG)
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(one * ufl.dx(domain=sub))),
                                         op=MPI.SUM)))
    t_cmp = time.perf_counter() - t0

    # 4. YAML side: the step file under its own scratch (gitignored), then the step block.
    t0 = time.perf_counter()
    src = _scratch("step-src")
    path = src / STEP_FILE
    if comm.rank == 0:
        m = _yaml_mapping()
        m["thermal"]["source_scale"] = 1.0 / float(ns["P"])
        m["thermal"]["transient"] = {"dt_s": DT_1 / 10.0, "n_steps": 1, "rho_c_j_per_m3_k": RHO_C,
                                     "theta": 1.0, "startup_half_steps": 0}
        path.write_text(yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    block = extract_tagged_block(text, STEP_TAG)
    assert block.count(f'"{STEP_FILE}"') == 1
    ns_y, out_y = _exec_block(monkeypatch, block.replace(f'"{STEP_FILE}"', repr(str(path)), 1),
                              "heating-step")
    t_file = time.perf_counter() - t0
    run = ns_y["run"]
    summ = run.summary["temperature"]

    # 4b. The page's CN study-file block, exactly as item 908 runs it
    # (`test_wf17_implant_rod_heating_tutorial.py:289–301`), under its own scratch.
    t0 = time.perf_counter()
    cn_src = _scratch("cn-file-src")
    cn_path = cn_src / "rod_heating.yaml"
    if comm.rank == 0:
        m = _yaml_mapping()
        m["thermal"]["source_scale"] = 1.0 / float(ns["P"])
        cn_path.write_text(yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    file_block = extract_tagged_block(text, HEAT_FILE_TAG)
    assert file_block.count('"rod_heating.yaml"') == 1
    ns_cn, out_cn = _exec_block(monkeypatch, file_block.replace('"rod_heating.yaml"',
                                                                repr(str(cn_path)), 1), "heating-file")
    t_cn = time.perf_counter() - t0
    run_cn = ns_cn["run"]
    summ_cn = run_cn.summary["temperature"]

    # 5. Rows.
    dt10 = DT_1 / 10.0
    th_s = {n: float(be10.balls[i]["theta_k"][1]) for i, n in enumerate(names)}
    th_y = {n: float(summ["balls"][i]["theta_k"][1]) for i, n in enumerate(names)}
    th_1 = {n: float(be.balls[i]["theta_k"][1]) for i, n in enumerate(names)}
    rows_s = _tip_slope_rows(th_s, dt10, comps, RHO_C)
    rows_y = _tip_slope_rows(th_y, dt10, comps, RHO_C)
    r_s, r_y, r_1 = _ratios(th_s, dt10, comps), _ratios(th_y, dt10, comps), _ratios(th_1, DT_1, comps)
    s_s = RHO_C * be10.volume_mean_k[1] / (dt10 * be10.source_power_w / v_sub)
    s_y = RHO_C * summ["volume_mean_k"][1] / (dt10 * summ["source_power_w"] / v_sub)

    th_cn, th_heat = _tips(summ_cn["balls"]), _tips(heat.balls)
    rel_theta = {n: abs(th_cn[n] / th_heat[n] - 1.0) for n in names}
    c_theta = abs(th_cn["top"] / th_heat["bottom"] - 1.0)
    vm_heat, vm_cn = float(heat.volume_mean_k[-1]), float(summ_cn["volume_mean_k"][-1])

    scale_ok = (summ["source_scale"] == cooling.source_scale and summ["n_steps"] == 1
                and summ["theta"] == 1.0)
    scale_cn_ok = summ_cn["source_scale"] == cooling.source_scale
    gates = {
        "Scale": scale_ok,
        "Scale CN": scale_cn_ok,
        "Cells": int(run.summary["cells"]) == ROD_EXACT_CELL_RECORD,
        "Cells CN": int(run_cn.summary["cells"]) == ROD_EXACT_CELL_RECORD,
        "E study": be10.max_balance_residual <= E_BAND,
        "E yaml": summ["max_balance_residual"] <= E_BAND,
        "E heat": heat.max_balance_residual <= E_BAND,
        "E CN": summ_cn["max_balance_residual"] <= E_BAND,
        "Band S2 study": abs(s_s - 1.0) <= S2_BAND,
        "Band S2 yaml": abs(s_y - 1.0) <= S2_BAND,
        "c-theta": c_theta > IDENTITY_BAND,
    }
    for n in names:
        gates[f"Y-theta {n}"] = rel_theta[n] <= IDENTITY_BAND
        gates[f"Band R_Pi study {n}"] = rows_s[n]["R_Pi"] <= S2_BAND
        gates[f"Band R_Pi yaml {n}"] = rows_y[n]["R_Pi"] <= S2_BAND
        gates[f"c-dt {n}"] = abs(r_y[n] - r_1[n]) > IDENTITY_BAND
        gates[f"R_Pi form {n}"] = abs(abs(r_y[n] - 1.0) - rows_y[n]["R_Pi"]) <= IDENTITY_BAND

    _p0(f"\n[WF-17 s43] -n {comm.size}: heating block {t_py:.1f} s, two BE steps {t_be:.1f} s, "
        f"comparand {t_cmp:.1f} s, step file block {t_file:.1f} s, CN file block {t_cn:.1f} s")
    _p0(f"[WF-17 s43] step block stdout (rank 0): {out_y!r}")
    _p0(f"[WF-17 s43] CN file block stdout (rank 0): {out_cn!r}")
    _p0(f"[WF-17 s43] (Scale) file source_scale {summ['source_scale']!r} vs block "
        f"{cooling.source_scale!r}; n_steps {summ['n_steps']}, theta {summ['theta']!r}, dt_s "
        f"{summ['dt_s']!r}: {'green' if scale_ok else 'RED'} (precondition, ASSERTED)")
    _p0(f"[WF-17 s43] (Scale CN) CN file source_scale {summ_cn['source_scale']!r} vs block "
        f"{cooling.source_scale!r}: {'green' if scale_cn_ok else 'RED'} (precondition, ASSERTED)")
    _p0(f"[WF-17 s43] (Cells) step file {run.summary['cells']}, CN file {run_cn.summary['cells']}, "
        f"block {heat.cells}, record {ROD_EXACT_CELL_RECORD} (ASSERTED)")
    _p0(f"[WF-17 s43] (E) CN heat block {heat.max_balance_residual:.3e}, CN file "
        f"{summ_cn['max_balance_residual']:.3e} (E_BAND {E_BAND:.0e}, ASSERTED)")
    for n in names:
        _p0(f"[WF-17 s43] (Y-theta) {n}: Study theta(t*) {th_heat[n]:.15e} K vs CN file "
            f"{th_cn[n]:.15e} K: rel {rel_theta[n]:.3e} vs item 908 {PRIOR_Y_THETA[n]} "
            f"(IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED)")
    _p0(f"[WF-17 s43] (c-theta) CN file top {th_cn['top']:.15e} K vs Study bottom "
        f"{th_heat['bottom']:.15e} K: rel {c_theta:.3e} > IDENTITY_BAND (miss ASSERTED; "
        f"predicted ~2.04e-02)")
    _p0(f"[WF-17 s43] volume-mean theta(t*) Study {vm_heat:.15e} K, CN file {vm_cn:.15e} K: rel "
        f"{abs(vm_cn / vm_heat - 1.0):.3e} (printed, no gate)")
    _p0(f"[WF-17 s43] P {float(np.real(ns['P'])):.12e} W vs item 884 {PRIOR_P} W (printed); "
        f"source_power_w Study {be10.source_power_w:.15e}, file {summ['source_power_w']:.15e}")
    _p0(f"[WF-17 s43] (E) Study BE dt1/10 {be10.max_balance_residual:.3e}, file "
        f"{summ['max_balance_residual']:.3e}, Study BE dt1 {be.max_balance_residual:.3e} "
        f"(E_BAND {E_BAND:.0e}, ASSERTED on the first two)")
    for n in names:
        _p0(f"[WF-17 s43] (Y-Pi) {n}: r_Y {r_y[n]!r} vs r_S {r_s[n]!r}: |dr| "
            f"{abs(r_y[n] - r_s[n]):.3e} beside (E) Study {be10.max_balance_residual:.3e} / file "
            f"{summ['max_balance_residual']:.3e} (reading, no band)")
        _p0(f"[WF-17 s43] (Band) {n}: R_Pi Study {rows_s[n]['R_Pi']:.4e}, file "
            f"{rows_y[n]['R_Pi']:.4e} vs item 884 {PRIOR_RPI_10[n]} (S2_BAND {S2_BAND:.0e}, ASSERTED)")
        _p0(f"[WF-17 s43] (c-dt) {n}: r_Y vs Study dt1 r {r_1[n]!r}: |dr| "
            f"{abs(r_y[n] - r_1[n]):.3e} > IDENTITY_BAND (miss ASSERTED); R_Pi at dt1 "
            f"{abs(r_1[n] - 1.0):.4e}")
        _p0(f"[WF-17 s43] block theta(t*) {n} {float(heat.balls[names.index(n)]['theta_k'][-1]):.12e} K "
            f"(item 908's reading, not re-gated)")
    _p0(f"[WF-17 s43] (Y-S2) s_Y {s_y!r} vs s_S {s_s!r}: |ds| {abs(s_y - s_s):.3e} beside (E) "
        f"Study {be10.max_balance_residual:.3e} / file {summ['max_balance_residual']:.3e} "
        f"(reading, no band)")
    _p0(f"[WF-17 s43] (Band) |s-1| Study {abs(s_s - 1.0):.4e}, file {abs(s_y - 1.0):.4e} vs item 884 "
        f"{PRIOR_S2_10} (S2_BAND {S2_BAND:.0e}, ASSERTED); V_sub {v_sub:.12e}")
    red = sorted(set(comm.allreduce([k for k, v in gates.items() if not v], op=MPI.SUM)))
    _p0(f"[WF-17 s43] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red
