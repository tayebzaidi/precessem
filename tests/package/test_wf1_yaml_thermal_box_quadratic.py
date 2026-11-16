"""`WF-1` step 16: ``excitation.source: {kind: transverse_quadratic}`` on a
``thermal_box`` file reproduces `MAT-21` step 9's (P36) records through
``femem run``.

Brief: docs/planning/chunks/WF-1.md, "§9 item 528" and the review ruling
2026-10-18 03:00 (point 7).  The solve legs are gated in the complex build (the
build of P36's records); they run only with ``FEM_EM_WF1_S16=1`` (and then
assert the complex build and ``-n 2``, never skip).  The parse node builds no
mesh and runs in either build.

Anchors: (K) the YAML Study = the Python Study (``_spec_diff``), cooling /
transient by ``float.hex``, ``femem check`` names the source class; (X)
refusals naming their field; (Pdeg) 116 025 P2 dofs on every solve; (I-Py)
the file's theta = the Python Study's at ``L0_BAND``; (S1-quad) / (T-quad)
interior at ``INTERIOR_BAND`` and at P36's ``.6e`` digits; (S2-quad) Gauss
at ``S2_GAUSS_BAND``; (E-quad) CN balance at ``E_BAND``; (U) the uniform P2
file reproduces step 3's m = 4 digits.  Controls (rule (e)): (c-src') and
(c-src) asserted, backed by P36's readings on the same fixture; (c-z)
printed only.
"""

from __future__ import annotations

import copy
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver import (
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
    UniformHeatSource,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study_yaml import _load, _study_from_mapping

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_study_yaml import LOGS, WF1, _cmp, _fresh_dir, _mapping, _rc_all, _spec_diff
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND
from tests.validation.test_mat21_box_closed_form import HALF, SERIES_TERMS, Q
from tests.validation.test_mat21_box_fem import T_STAR, _interior_chi, _steady_ex, _trans_ex
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _exact_p, _rel_p
from tests.validation.test_mat21_box_quadratic_closed_form import WEIGHTS_NOZ
from tests.validation.test_mat21_box_quadratic_fem import _quad_ex, _quad_trans_ex
from tests.validation.test_mat21_box_separable_fem import STEP3_M4

S16_ENV = "FEM_EM_WF1_S16"
QUAD_STEADY = WF1 / "thermal_box_quadratic_steady_m4.yaml"
QUAD_CN = WF1 / "thermal_box_quadratic_cn_m4.yaml"
UNIFORM_P2 = WF1 / "thermal_box_uniform_steady_m4_p2.yaml"
BALL = WF1 / "thermal_ball_duty_a.yaml"
DIV4 = (52, 32, 8)

# P36's records, transcribed from 20261016T171904Z_MAT-21-step9-W2-quadratic-fem.log.
P36_S1 = "7.219878e-05"  # :76 (S1-quad) steady interior
P36_T = "2.694340e-04"  # :83 (T-quad) CN interior at t*
P36_E = 2.285e-13  # :84 (E-quad), printed beside
P36_S2 = 9.775e-14  # :73 (S2-quad), printed beside
P36_DOFS = 116025  # :70 (Pdeg)
P36_C_SRCP = 5.1471e-01  # :78 (c-src')
P36_C_SRC = 4.4856e-01  # :77 (c-src)
P36_C_Z = 4.7893e-02  # :79 (c-z), printed only


def _knob() -> bool:
    raw = os.environ.get(S16_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S16_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _all(comm, flag):
    return bool(comm.allreduce(bool(flag), op=MPI.LAND))


def _field(fn):
    try:
        fn()
    except StudyConfigError as exc:
        return exc.field
    return None


def _with(m, value, *path):
    m = copy.deepcopy(m)
    d = m
    for p in path[:-1]:
        d = d[p]
    d[path[-1]] = value
    return m


def _quad_study():
    return Study(ThermalBoxGeometry(HALF, DIV4), None, TransverseQuadraticHeatSource(Q))


def _cooling():
    return RobinCooling(K_GEL, H_GEL, T_INF, degree=DEGREE)


# ------------------------------------------------------------------ no mesh


def test_wf1_s16_parse_and_refusals(monkeypatch, capsys):
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

    # (K) the Study spec
    want_spec = {
        QUAD_STEADY: _quad_study(),
        QUAD_CN: _quad_study(),
        UNIFORM_P2: Study(ThermalBoxGeometry(HALF, DIV4), None, UniformHeatSource(Q)),
    }
    diffs = {f.name: _spec_diff(load_study(f), s) for f, s in want_spec.items()}
    _say(capsys, f"\n[WF-1 s16] (K) spec diffs {diffs}")
    assert all(d == [] for d in diffs.values()), diffs
    assert type(load_study(UNIFORM_P2).excitation) is UniformHeatSource
    assert type(load_study(QUAD_STEADY).excitation) is TransverseQuadraticHeatSource

    # (K) the fourth element by hex: cooling at degree 2; P36's CN dt = T_STAR / N_MAIN
    want = {
        QUAD_STEADY: (_cooling(), None),
        UNIFORM_P2: (_cooling(), None),
        QUAD_CN: (_cooling(), TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)),
    }
    for f, (cool, tr) in want.items():
        _, kind, outputs, th = _load(f)
        out = []
        _cmp(th.cooling, cool, "cooling", out)
        _cmp(th.transient, tr, "transient", out)
        _say(capsys, f"  (K) {f.name}: kind {kind!r}, outputs {outputs!r}, thermal diffs {out}")
        assert kind == "thermal_box"
        assert th.duty is None
        assert out == [], (f.name, out)
    th_c = _load(QUAD_CN)[3]
    _say(capsys, f"  (K) CN dt_s {th_c.transient.dt_s.hex()} (T_STAR / N_MAIN "
         f"{(T_STAR / N_MAIN).hex()})")

    # (K) femem check names the source class
    for f, name in ((QUAD_STEADY, "TransverseQuadraticHeatSource"),
                    (QUAD_CN, "TransverseQuadraticHeatSource"),
                    (UNIFORM_P2, "UniformHeatSource")):
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        lines = capsys.readouterr().out.splitlines()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}; rank-0 lines {lines}")
        assert rc_all == (0, 0)
        if comm.rank == 0:
            assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]
            assert f"excitation: {name}" in lines

    # (X) the variant table
    s = _mapping(QUAD_STEADY)
    rows = [
        (_with(s, {"kind": "parabolic"}, "excitation", "source"), "excitation.source.kind"),
        (_with(s, [1], "excitation", "source"), "excitation.source"),
        (_with(s, {"kind": "transverse_quadratic", "weights": [1, 2]}, "excitation", "source"),
         "excitation.source.weights"),
        (_with(s, {}, "excitation", "source"), "excitation.source.kind"),
        (_with(_mapping(BALL), {"kind": "transverse_quadratic"}, "excitation", "source"),
         "excitation.source"),
    ]
    got = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for (_m, want_f), g in zip(rows, got):
        _say(capsys, f"  (X) want {want_f!r} got {g!r}")
    assert got == [f for _, f in rows], [(f, g) for (_, f), g in zip(rows, got) if f != g]
    # `source: {kind: uniform}` is the absent key's Study.
    d_u = _spec_diff(_study_from_mapping(_with(_mapping(UNIFORM_P2), {"kind": "uniform"},
                                                "excitation", "source")),
                     load_study(UNIFORM_P2))
    _say(capsys, f"  (K) explicit uniform vs absent source: spec diff {d_u}")
    assert d_u == []

    # The API: a ball refuses the new class; q must be positive.
    with pytest.raises(TypeError):
        Study(ThermalBallGeometry(0.05, 0.01), None, TransverseQuadraticHeatSource(Q))
    with pytest.raises(ValueError):
        TransverseQuadraticHeatSource(0.0)
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _solve_preamble():
    if not _knob():
        pytest.skip(f"{S16_ENV} unset: the quadratic-box solve legs are environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "P36's records are complex-build records"
    assert comm.size == 2, comm.size
    return comm


def _ndofs(T):
    return int(T.function_space.dofmap.index_map.size_global)


def _on_p(msh_ref, T):
    """``T`` (P2 on a bitwise-identical mesh build) as a Function on ``msh_ref``:
    the P2 analogue of ``test_wf15_study_transient._on``, which is P1-only."""
    comm = msh_ref.comm
    V = fem.functionspace(msh_ref, ("Lagrange", DEGREE))
    same = np.array_equal(
        V.tabulate_dof_coordinates(), T.function_space.tabulate_dof_coordinates()
    ) and V.dofmap.index_map.size_local == T.function_space.dofmap.index_map.size_local
    assert _all(comm, same), "P2 dof layout differs between the two builds"
    f = fem.Function(V)
    f.x.array[:] = T.x.array
    f.x.scatter_forward()
    return f


def _owned_rel(T_a, T_b, comm):
    """Owned max |T_a − T_b| / max |T_b − T∞|, reduced (bs 1)."""
    n = T_b.function_space.dofmap.index_map.size_local
    a = np.asarray(T_a.x.array[:n])
    b = np.asarray(T_b.x.array[:n])
    dmax = comm.allreduce(float(np.max(np.abs(a - b), initial=0.0)), op=MPI.MAX)
    nmax = comm.allreduce(float(np.max(np.abs(b - T_INF), initial=0.0)), op=MPI.MAX)
    return dmax / nmax


def test_wf1_s16_quadratic_steady_reproduces_p36():
    comm = _solve_preamble()
    t00 = time.perf_counter()
    out = LOGS / "thermal_box_quadratic_steady"
    _fresh_dir(out, comm)
    run = run_study_file(QUAD_STEADY, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    th = run.thermal
    msh = th.submesh
    py = _quad_study().temperature(_cooling())
    T_py = _on_p(msh, py.temperature)
    ip = _owned_rel(T_py, th.temperature, comm)
    out_u = LOGS / "thermal_box_uniform_p2"
    _fresh_dir(out_u, comm)
    run_u = run_study_file(UNIFORM_P2, output_dir=out_u, comm=comm)
    T_u = _on_p(msh, run_u.thermal.temperature)
    t_solves = time.perf_counter() - t00
    dofs = [_ndofs(th.temperature), _ndofs(py.temperature), _ndofs(run_u.thermal.temperature)]

    th_q = th.temperature - T_INF
    th_u = T_u - T_INF
    chi = _interior_chi(msh)
    # (S2-quad) P36's Gauss reading
    ds = ufl.Measure("ds", domain=msh, subdomain_data=th.facet_tags)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_q * ds)))
    sp = float(np.real(th.source_power))
    gauss = abs(conv - sp) / sp
    # (U) step 3's digits through the uniform P2 file
    ex_u, npts = _exact_p(msh, _steady_ex())
    u_whole = _rel_p(msh, th_u, ex_u)
    u_int = _rel_p(msh, th_u, ex_u, chi)
    # (S1-quad)
    ex_q, _ = _exact_p(msh, _quad_ex())
    s1 = _rel_p(msh, th_q, ex_q, chi)
    # controls
    c_srcp = _rel_p(msh, th_u, ex_q, chi)
    c_src = _rel_p(msh, th_q, ex_u, chi)
    ex_z, _ = _exact_p(msh, _quad_ex(WEIGHTS_NOZ, SERIES_TERMS))
    c_z = _rel_p(msh, th_q, ex_z, chi)

    _p0(f"WF1S16 (Pdeg) steady P2 dofs file / Python / uniform file {dofs} (P36 {P36_DOFS}); "
        f"file run {t_file:.1f}s, three solves {t_solves:.1f}s; summary source "
        f"{run.summary['temperature']['source']!r} / {run_u.summary['temperature']['source']!r}")
    _p0(f"WF1S16 (I-Py) steady: owned max|dθ|/max|θ| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e})")
    _p0(f"WF1S16 (S1-quad) steady interior vs quadratic series {s1:.6e} (P36 {P36_S1}; band "
        f"{INTERIOR_BAND:.0e}); {npts} P3 points")
    _p0(f"WF1S16 (S2-quad) Gauss |h oint theta - source_power|/source_power = {gauss:.3e} (P36 "
        f"{P36_S2:.3e}; band {S2_GAUSS_BAND:.0e}); source_power/(Q V) - 1 = "
        f"{sp / (Q * 8.0 * HALF[0] * HALF[1] * HALF[2]) - 1.0:+.3e}")
    _p0(f"WF1S16 (U) uniform P2 file: whole {u_whole:.6e} (step 3 {STEP3_M4['whole']}); interior "
        f"{u_int:.6e} (step 3 {STEP3_M4['interior']})")
    _p0(f"WF1S16 control (c-src') uniform file vs quadratic series {c_srcp:.4e} (P36 "
        f"{P36_C_SRCP:.4e}; asserted > {INTERIOR_BAND:.0e})")
    _p0(f"WF1S16 control (c-src) quadratic file vs uniform series {c_src:.4e} (P36 "
        f"{P36_C_SRC:.4e}; asserted > {INTERIOR_BAND:.0e})")
    _p0(f"WF1S16 control (c-z) quadratic file vs no-z series {c_z:.4e} (P36 {P36_C_Z:.4e}; "
        f"printed only, ratio to band {c_z / INTERIOR_BAND:.2f})")
    assert dofs == [P36_DOFS] * 3
    assert run.summary["temperature"]["source"] == "transverse_quadratic"
    assert run_u.summary["temperature"]["source"] == "uniform"
    assert run.summary["temperature"]["thermal_degree"] == DEGREE
    assert ip <= L0_BAND, f"(I-Py) steady {ip:.3e}"
    assert s1 <= INTERIOR_BAND
    assert f"{s1:.6e}" == P36_S1
    assert gauss <= S2_GAUSS_BAND
    assert f"{u_whole:.6e}" == STEP3_M4["whole"]
    assert f"{u_int:.6e}" == STEP3_M4["interior"]
    assert c_srcp > INTERIOR_BAND
    assert c_src > INTERIOR_BAND
    _p0(f"WF1S16 W1a leg {time.perf_counter() - t00:.1f}s")


def test_wf1_s16_quadratic_cn_reproduces_p36():
    comm = _solve_preamble()
    t00 = time.perf_counter()
    out = LOGS / "thermal_box_quadratic_cn"
    _fresh_dir(out, comm)
    run = run_study_file(QUAD_CN, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    th = run.thermal
    msh = th.mesh
    py = _quad_study().temperature(_cooling(),
                                   transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))
    T_py = _on_p(msh, py.temperature)
    ip = _owned_rel(T_py, th.temperature, comm)
    dofs = [_ndofs(th.temperature), _ndofs(py.temperature)]
    mult = np.array_equal(np.asarray(th.source_multipliers), np.asarray(py.source_multipliers))
    th_t = th.temperature - T_INF
    chi = _interior_chi(msh)
    ex_t, npts = _exact_p(msh, _quad_trans_ex(T_STAR))
    t_int = _rel_p(msh, th_t, ex_t, chi)
    ex_ut, _ = _exact_p(msh, _trans_ex(T_STAR))
    c_src_t = _rel_p(msh, th_t, ex_ut, chi)
    bal = run.summary["temperature"]["max_balance_residual"]
    _p0(f"WF1S16 (Pdeg) CN P2 dofs file / Python {dofs} (P36 {P36_DOFS}); multipliers equal "
        f"{mult}; file run {t_file:.1f}s; summary source {run.summary['temperature']['source']!r}")
    _p0(f"WF1S16 (I-Py) CN: owned max|dθ|/max|θ| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e})")
    _p0(f"WF1S16 (T-quad) CN N={N_MAIN} at t*: interior vs quadratic series {t_int:.6e} (P36 "
        f"{P36_T}; band {INTERIOR_BAND:.0e}); {npts} P3 points")
    _p0(f"WF1S16 (E-quad) CN max balance residual {bal:.3e} (P36 {P36_E:.3e}; band {E_BAND:.0e})")
    _p0(f"WF1S16 control (c-src-T) CN quadratic vs uniform transient series {c_src_t:.4e} "
        f"(P36 6.9562e-01; printed only)")
    assert dofs == [P36_DOFS] * 2
    assert mult
    assert run.summary["temperature"]["source"] == "transverse_quadratic"
    assert run.summary["temperature"]["duty"] is None
    assert ip <= L0_BAND, f"(I-Py) CN {ip:.3e}"
    assert t_int <= INTERIOR_BAND
    assert f"{t_int:.6e}" == P36_T
    assert bal <= E_BAND
    _p0(f"WF1S16 W1b leg {time.perf_counter() - t00:.1f}s")
