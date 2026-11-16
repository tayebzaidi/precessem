"""`WF-1` step 18: a ``thermal_box`` file with ``excitation.source:
{kind: transverse_quadratic}`` and ``thermal.duty`` reproduces `MAT-21`
step 11's (P38) duty-cycled transverse-source record through ``femem run``.

Brief: docs/planning/chunks/WF-1.md, "§9 item 540".  The solve leg is gated in
the complex build (the build of P38's record); it runs only with
``FEM_EM_WF1_S18=1`` (and then asserts the complex build and ``-n 2``, never
skips).  The parse node builds no mesh and runs in either build.

Anchors: (K) the YAML Study = the Python Study (``_spec_diff``), cooling /
transient / duty by ``float.hex``, ``femem check`` names the source class;
(X) the file without ``transient`` raises at ``thermal.duty``; (Pdeg) 116 025
P2 dofs on both solves; (UT) the file's mesh = the Python Study's; (S0) the
schedule string; (I-Py) theta_file = theta_Python at ``L0_BAND``; (T-qduty)
interior at ``INTERIOR_BAND`` and at P38's ``.6e`` digits; (E-qduty) CN
balance at ``E_BAND``; (I=) the summary is its own result.  Control (rule
(e)): (c-step) asserted, backed by P38's reading on the same fixture.
"""

from __future__ import annotations

import copy
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    DutyCycle,
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    StudyConfigError,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    TransverseQuadraticHeatSource,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study_yaml import _load, _study_from_mapping

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_study_yaml import LOGS, WF1, _cmp, _fresh_dir, _mapping, _rc_all, _spec_diff
from tests.package.test_wf1_yaml_thermal_box import _check_summary
from tests.package.test_wf1_yaml_thermal_box_quadratic import P36_DOFS, _on_p, _owned_rel
from tests.package.test_wf15_study_transient import _same_mesh
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF
from tests.validation.test_mat21_box_closed_form import HALF, Q
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, DUTY, HALF_T, N_STEPS
from tests.validation.test_mat21_box_fem import T_STAR, _interior_chi
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _exact_p, _rel_p
from tests.validation.test_mat21_box_quadratic_duty_closed_form import _quad_duty
from tests.validation.test_mat21_box_quadratic_duty_fem import SCHEDULE
from tests.validation.test_mat21_box_quadratic_fem import _quad_trans_ex

S18_ENV = "FEM_EM_WF1_S18"
QUAD_DUTY = WF1 / "thermal_box_quadratic_duty_m4.yaml"
DIV4 = (52, 32, 8)
CELLS_M4 = 79872  # (UT) the m = 4 box, as the brief records

# P38's records, transcribed from 20261018T124640Z_MAT-21-step11-W2-quadratic-duty-fem.log.
P38_T_QDUTY = "2.543310e-04"  # :82 (T-qduty) interior at t*
P38_E = 2.288e-13  # :83 (E-qduty), printed beside
P38_C_STEP = 5.0333e-01  # :84 (c-step)
P38_DOFS = 116025  # :89 (Pdeg)


def _knob() -> bool:
    raw = os.environ.get(S18_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S18_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _field(fn):
    try:
        fn()
    except StudyConfigError as exc:
        return exc.field
    return None


def _without(m, *path):
    m = copy.deepcopy(m)
    d = m
    for p in path[:-1]:
        d = d[p]
    del d[path[-1]]
    return m


def _quad_study():
    return Study(ThermalBoxGeometry(HALF, DIV4), None, TransverseQuadraticHeatSource(Q))


def _cooling():
    return RobinCooling(K_GEL, H_GEL, T_INF, degree=DEGREE)


def _transient():
    return TransientSpec(DT_BOX, N_STEPS, RHO_C)


def _duty():
    return DutyCycle(2.0 * HALF_T, DUTY)


# ------------------------------------------------------------------ no mesh


def test_wf1_s18_parse_and_refusals(monkeypatch, capsys):
    """(K), (X): no mesh, no solve, either build."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                ThermalBallGeometry, ThermalBoxGeometry):
        monkeypatch.setattr(cls, "build", stub)

    assert P38_DOFS == P36_DOFS  # same fixture (m = 4, P2)

    # (K) the Study spec
    d = _spec_diff(load_study(QUAD_DUTY), _quad_study())
    _say(capsys, f"\n[WF-1 s18] (K) spec diff {d}")
    assert d == []
    assert type(load_study(QUAD_DUTY).excitation) is TransverseQuadraticHeatSource

    # (K) cooling / transient / duty by hex
    _, kind, outputs, th = _load(QUAD_DUTY)
    out = []
    _cmp(th.cooling, _cooling(), "cooling", out)
    _cmp(th.transient, _transient(), "transient", out)
    _cmp(th.duty, _duty(), "duty", out)
    _say(capsys, f"  (K) {QUAD_DUTY.name}: kind {kind!r}, outputs {outputs!r}, thermal diffs {out}")
    _say(capsys, f"  (K) dt_s {th.transient.dt_s.hex()} (DT_BOX {DT_BOX.hex()}); period "
         f"{th.duty.period_s.hex()} (2 HALF_T {(2.0 * HALF_T).hex()}); duty "
         f"{th.duty.duty.hex()} (DUTY {float(DUTY).hex()})")
    assert kind == "thermal_box"
    assert out == [], out

    # (K) femem check names the source class
    capsys.readouterr()
    rc = cli.main(["check", str(QUAD_DUTY)])
    lines = capsys.readouterr().out.splitlines()
    rc_all = _rc_all(rc, comm)
    _say(capsys, f"  (K) femem check {QUAD_DUTY.name}: rc {rc_all}; rank-0 lines {lines}")
    assert rc_all == (0, 0)
    if comm.rank == 0:
        assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]
        assert "excitation: TransverseQuadraticHeatSource" in lines

    # (X) the same file without `transient` raises at thermal.duty
    got = _field(lambda: _study_from_mapping(_without(_mapping(QUAD_DUTY), "thermal", "transient")))
    _say(capsys, f"  (X) quadratic duty file without transient: want 'thermal.duty' got {got!r}")
    assert got == "thermal.duty"
    assert calls["n"] == 0


# ------------------------------------------------------------------ solve


def _ndofs(T):
    return int(T.function_space.dofmap.index_map.size_global)


def _hex_dict_equal(a, b):
    return set(a) == set(b) and all(float(a[k]).hex() == float(b[k]).hex() for k in a)


def test_wf1_s18_quadratic_duty_reproduces_p38():
    if not _knob():
        pytest.skip(f"{S18_ENV} unset: the quadratic-duty solve leg is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "P38's records are complex-build records"
    assert comm.size == 2, comm.size

    t00 = time.perf_counter()
    out = LOGS / "thermal_box_quadratic_duty"
    _fresh_dir(out, comm)
    run = run_study_file(QUAD_DUTY, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    th = run.thermal
    msh = th.mesh
    py = _quad_study().temperature(_cooling(), transient=_transient(), duty=_duty())
    t_py = time.perf_counter() - t00 - t_file
    ca, cb, geo = _same_mesh(msh, py.mesh, comm)
    T_py = _on_p(msh, py.temperature)
    ip = _owned_rel(T_py, th.temperature, comm)
    dofs = [_ndofs(th.temperature), _ndofs(py.temperature)]
    st = run.summary["temperature"]
    s0 = "".join(str(int(m)) for m in st["source_multipliers"])
    mult = np.array_equal(np.asarray(th.source_multipliers), np.asarray(py.source_multipliers))

    th_t = th.temperature - T_INF
    chi = _interior_chi(msh)
    ex_d, npts = _exact_p(msh, lambda x, y, z: _quad_duty(x, y, z, T_STAR))
    t_qduty = _rel_p(msh, th_t, ex_d, chi)
    ex_s, _ = _exact_p(msh, _quad_trans_ex(T_STAR))
    c_step = _rel_p(msh, th_t, ex_s, chi)
    bal = st["max_balance_residual"]

    _p0(f"WF1S18 (Pdeg) P2 dofs file / Python {dofs} (P38 {P38_DOFS}); file run {t_file:.1f}s, "
        f"Python run {t_py:.1f}s; summary source {st['source']!r}; duty {st['duty']!r}")
    _p0(f"WF1S18 (UT) cells file / Python {ca} / {cb} (record {CELLS_M4}); geometry bitwise {geo}")
    _p0(f"WF1S18 (S0) schedule {s0} (SCHEDULE {SCHEDULE}); multipliers equal Python {mult}")
    _p0(f"WF1S18 (I-Py) owned max|dθ|/max|θ| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e})")
    _p0(f"WF1S18 (T-qduty) interior vs duty quadratic series {t_qduty:.6e} (P38 {P38_T_QDUTY}; "
        f"band {INTERIOR_BAND:.0e}); {npts} P3 points")
    _p0(f"WF1S18 (E-qduty) CN max balance residual {bal:.3e} (P38 {P38_E:.3e}; band {E_BAND:.0e})")
    _p0(f"WF1S18 control (c-step) duty file vs unscheduled transient series {c_step:.4e} "
        f"(P38 {P38_C_STEP:.4e}; asserted > {INTERIOR_BAND:.0e})")
    _check_summary(run, comm, "quadratic duty", transient=True)

    assert dofs == [P38_DOFS] * 2
    assert ca == cb == CELLS_M4 and geo
    assert s0 == SCHEDULE
    assert mult
    assert st["source"] == "transverse_quadratic"
    assert _hex_dict_equal(st["duty"], {"period_s": 2.0 * HALF_T, "duty": DUTY})
    assert ip <= L0_BAND, f"(I-Py) {ip:.3e}"
    assert t_qduty <= INTERIOR_BAND
    assert f"{t_qduty:.6e}" == P38_T_QDUTY
    assert bal <= E_BAND
    assert c_step > INTERIOR_BAND
    _p0(f"WF1S18 W1 leg {time.perf_counter() - t00:.1f}s")
