"""`WF-15` step 3 (§9 item 330, chain P16): ``Study.temperature`` with a duty
cycle and a container wall, on thermal-only balls (no EM solve).

Brief: docs/planning/chunks/WF-15.md, "§9 item 330" and its review rulings
(2026-10-08 03:00).  Real build (the records' build).  Records reproduced:
`MAT-17` step 2 cases A / B (``20261003T171113Z_MAT-17-step2-W1b.log:98, 101``)
and `MAT-18` step 2(b) (T1) (``20261003T184850Z_MAT-18-step2b-W1b.log:88``),
citations verified against the logs 2026-10-08 before registering.

Anchors: (UT) bitwise upstream of the solve (cells, mesh geometry, DG0 k / ρc
/ q, the schedule's multipliers) with (S0); (IT) API vs the direct
``solve_transient_conduction`` ≤ ``L0_BAND``; (G) interior ≤ 1 % and balance
≤ 1e-9 (`MAT-9`'s literals, by value); (R) the records at ``.4e``.
Controls (rule (e)): (c1) duty dropped, asserted > 1 %, factor **predicted**
≈ 1.0 (ruling (6)); (c2) wall ρc dropped, asserted > 1 %, prior
same-fixture measurement 5.3140e-02 (``…step2b-W1b.log:94``); (c3) planted
``np.spacing`` read back exactly (asserted, deterministic).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver import (
    ContainerWall,
    DutyCycle,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    TransientSpec,
    UniformHeatSource,
)
from fem_em_solver.utils.analytical import robin_sphere_duty_cycle_temperature
from tests.materials.test_phantom_material_model import _dg0_owned_diff
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat9_transient_conduction_sphere import (
    RHO_C,
    T_INF,
    H,
    K,
    Q,
    R,
    _ball,
    _interior_error,
)
from tests.validation.test_mat17_duty_cycle_closed_form import CASES
from tests.validation.test_mat17_duty_cycle_fem import (
    DT_B,
    ON_TIME_A,
    ON_TIME_B,
    T_END_A,
    T_END_B,
    _census,
    _chi,
    _n,
)
from tests.validation.test_mat17_duty_cycle_fem import _run as _run17
from tests.validation.test_mat18_two_layer_sphere_closed_form import (
    A_,
    B_,
    K1,
    K2,
    Q1,
    RC1,
    RC2,
)
from tests.validation.test_mat18_two_layer_sphere_transient_fem import (
    CELL_GUARD,
    N_CN,
    T_STAR,
    _errors,
    _fine,
)
from tests.validation.test_mat18_two_layer_sphere_transient_fem import _run as _run18

P_A, D_A = CASES["A"]
P_B, D_B = CASES["B"]
RES = 0.005
CHUNK = 20_000

# `MAT-9` step 2's bands, by value (test_mat9_transient_conduction_sphere.py
# :274–275, 301–302, 321–322).
INTERIOR_BAND = 1e-2
BALANCE_BAND = 1e-9

# Records (verified 2026-10-08 against the cited lines).
MAT17_A_PRINT = "4.3870e-04"  # 20261003T171113Z_MAT-17-step2-W1b.log:98
MAT17_B_PRINT = "1.0683e-03"  # 20261003T171113Z_MAT-17-step2-W1b.log:101
MAT18_T1_PRINT = "8.2066e-04"  # 20261003T184850Z_MAT-18-step2b-W1b.log:88
MAT18_T1_CELLS = 215542  # 20261003T184850Z_MAT-18-step2b-W1b.log:88
MAT18_C1_UNIFORM_RHOC = 5.3140e-02  # 20261003T184850Z_MAT-18-step2b-W1b.log:94
C1_PREDICTED = 1.0  # ruling (6): predicted, no prior same-fixture measurement


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _rss_gib(comm):
    import resource

    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(kib, op=MPI.SUM) / 1024.0**2


def _all(comm, flag):
    return bool(comm.allreduce(bool(flag), op=MPI.LAND))


def _same_mesh(m_a, m_b, comm):
    """Bitwise: global cells, and the rank-local geometry arrays (reduced)."""
    ca = m_a.topology.index_map(3).size_global
    cb = m_b.topology.index_map(3).size_global
    geo = np.array_equal(m_a.geometry.x, m_b.geometry.x) and np.array_equal(
        m_a.geometry.dofmap, m_b.geometry.dofmap
    )
    return ca, cb, _all(comm, geo)


def _on(msh_ref, T):
    """``T`` (P1 on a bitwise-identical mesh build) as a Function on ``msh_ref``.

    UFL forms cannot mix two mesh objects; the copy is exact because the two
    builds' geometry, and so the P1 dof layout and coordinates, are asserted
    identical first."""
    comm = msh_ref.comm
    V = fem.functionspace(msh_ref, ("Lagrange", 1))
    same = (
        np.array_equal(V.tabulate_dof_coordinates(), T.function_space.tabulate_dof_coordinates())
        and V.dofmap.index_map.size_local == T.function_space.dofmap.index_map.size_local
    )
    assert _all(comm, same), "P1 dof layout differs between the two builds"
    f = fem.Function(V)
    f.x.array[:] = T.x.array
    f.x.scatter_forward()
    return f


def _theta_owned_rel(T_a, T_b, comm):
    """Owned max |T_a − T_b| / max |T_b − T∞|, reduced (P1, bs 1)."""
    n = T_b.function_space.dofmap.index_map.size_local
    a = np.asarray(T_a.x.array[:n])
    b = np.asarray(T_b.x.array[:n])
    dmax = comm.allreduce(float(np.max(np.abs(a - b), initial=0.0)), op=MPI.MAX)
    nmax = comm.allreduce(float(np.max(np.abs(b - T_INF), initial=0.0)), op=MPI.MAX)
    return dmax / nmax


def _exact_duty(msh, t, period, duty):
    """θ(t) of `MAT-17` step 1's closed form at the P2 dof coordinates
    (owned + ghost), numpy, chunks of 20 000 distinct radii (brief trap 1)."""
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    f = fem.Function(V2)
    x = V2.tabulate_dof_coordinates()
    r_all = np.sqrt(np.sum(x**2, axis=1))
    r, inv = np.unique(r_all, return_inverse=True)
    out = np.empty(r.shape)
    for a in range(0, r.size, CHUNK):
        rr = r[a : a + CHUNK]
        out[a : a + CHUNK] = (
            robin_sphere_duty_cycle_temperature(
                rr,
                np.full(rr.shape, t),
                radius=R,
                k=K,
                rho_c=RHO_C,
                q=Q,
                h=H,
                t_inf=T_INF,
                period=period,
                duty=duty,
                periodic=False,
            )
            - T_INF
        )
    f.x.array[:] = out[inv.reshape(-1)]
    f.x.scatter_forward()
    return f


def _duty_study(dt, n, duty):
    st = Study(ThermalBallGeometry(R, RES), None, UniformHeatSource(Q))
    return st.temperature(
        RobinCooling(K, H, T_INF), transient=TransientSpec(dt, n, RHO_C), duty=duty
    )


# ------------------------------------------------------------------ no solve


def test_wf15_transient_spec_and_refusals():
    spec = TransientSpec(150.0, 94, RHO_C)
    assert (spec.theta, spec.startup_half_steps, spec.t_init_c) == (0.5, 2, None)
    geo = ThermalBallGeometry(A_, RES, wall=ContainerWall(B_ - A_, K2, RC2))
    assert geo.outer_radius == B_
    assert ThermalBallGeometry(R, RES).outer_radius == R
    # The schedule is `MAT-17`'s χ, value for value, on a grid that hits every
    # switch time of case A and case B.
    for P, D in (CASES["A"], CASES["B"]):
        d, chi = DutyCycle(P, D), _chi(P, D)
        ts = np.concatenate(
            [np.linspace(0.0, 13 * P, 4001), np.arange(14) * P, np.arange(14) * P + D * P]
        )
        assert all(d(float(t)) == chi(float(t)) for t in ts)
    for bad in (
        lambda: TransientSpec(0.0, 10, RHO_C),
        lambda: TransientSpec(1.0, 0, RHO_C),
        lambda: TransientSpec(1.0, 10, RHO_C, theta=1.5),
        lambda: TransientSpec(1.0, 10, RHO_C, startup_half_steps=1),
        lambda: DutyCycle(100.0, 0.0),
        lambda: DutyCycle(100.0, 1.5),
        lambda: DutyCycle(-1.0, 0.5),
        lambda: ContainerWall(0.0, K2, RC2),
        lambda: UniformHeatSource(0.0),
        lambda: ThermalBallGeometry(0.0, RES),
    ):
        with pytest.raises(ValueError):
            bad()
    with pytest.raises(TypeError):
        Study(ThermalBallGeometry(R, RES), "gel", UniformHeatSource(Q))
    with pytest.raises(TypeError):
        Study(ThermalBallGeometry(R, RES), None, None)
    st = Study(ThermalBallGeometry(R, RES), None, UniformHeatSource(Q))
    with pytest.raises(TypeError):
        st.solve()
    with pytest.raises(TypeError):
        st.temperature(RobinCooling(K, H, T_INF))  # no TransientSpec
    with pytest.raises(TypeError):
        st.temperature(RobinCooling(K, H, T_INF), transient=spec, duty=(P_A, D_A))
    with pytest.raises(ValueError):
        st.temperature(RobinCooling(K, H, T_INF, source_scale=2.0), transient=spec)
    assert st.result is None
    _say("WF15S spec defaults and refusals: green (no solve)")


# --------------------------------------------------------------- W1a, duty


def test_wf15_duty_cycle_reproduces_mat17():
    comm = MPI.COMM_WORLD
    t00 = time.perf_counter()
    msh, _, ft = _ball(RES)
    n_glob = msh.topology.index_map(3).size_global
    rows = {}
    for lab, P, D, dt, t_end, on_time, rec in (
        ("A", P_A, D_A, 150.0, T_END_A, ON_TIME_A, MAT17_A_PRINT),
        ("B", P_B, D_B, DT_B, T_END_B, ON_TIME_B, MAT17_B_PRINT),
    ):
        n = _n(t_end, dt)
        t0 = time.perf_counter()
        api = _duty_study(dt, n, DutyCycle(P, D))
        t_api = time.perf_counter() - t0
        direct = _run17(msh, ft, dt, n, 0.5, 2, schedule=_chi(P, D))
        # (UT) + (S0)
        ca, cb, geo = _same_mesh(api.mesh, msh, comm)
        mult_eq = np.array_equal(
            np.asarray(api.source_multipliers), np.asarray(direct.source_multipliers)
        )
        _census(api, dt, 2, on_time, f"WF15S (S0) case {lab} API")
        T_api = _on(msh, api.temperature)
        it = _theta_owned_rel(T_api, direct.temperature, comm)
        ex = _exact_duty(msh, t_end, P, D)
        e = _interior_error(msh, T_api, ex)
        e_dir = _interior_error(msh, direct.temperature, ex)
        rows[lab] = (e, ex, n, dt)
        _say(
            f"WF15S (UT) case {lab}: cells {ca} / {cb} (direct {n_glob}); geometry equal {geo}; "
            f"multipliers equal {mult_eq} ({len(api.source_multipliers)} (half-)steps)"
        )
        _say(
            f"WF15S (IT) case {lab}: owned max|dθ|/max|θ| = {it:.3e} (band {L0_BAND:.0e}, "
            f"predicted <= 1e-14)"
        )
        _say(
            f"WF15S (G) case {lab} CN+Rannacher dt={dt:.4f} N={n}: interior_rel_L2={e:.10e} "
            f"(band 1e-2) direct={e_dir:.10e}; API max_balance={api.max_balance_residual:.3e} "
            f"direct={direct.max_balance_residual:.3e}; API solve+mesh {t_api:.1f}s"
        )
        _say(f"WF15S (R) case {lab}: {e:.4e} vs record {rec}")
        assert ca == cb == n_glob, f"(UT) cells case {lab}"
        assert geo, f"(UT) geometry case {lab}"
        assert mult_eq, f"(UT) multipliers case {lab}"
        assert it <= L0_BAND, f"(IT) case {lab} {it:.3e}"
        assert api.max_balance_residual <= BALANCE_BAND, f"(G) balance case {lab}"
        assert e <= INTERIOR_BAND, f"(G) case {lab} {e:.4e}"
        assert f"{e:.4e}" == rec, f"(R) case {lab} {e:.4e} != {rec}"

    # (c1) the schedule dropped through the API (no DutyCycle), case A.
    e_a, ex_a, n_a, dt_a = rows["A"]
    c1 = _duty_study(dt_a, n_a, None)
    e_c1 = _interior_error(msh, _on(msh, c1.temperature), ex_a)
    _say(
        f"WF15S (c1) no DutyCycle, case A: interior miss measured {e_c1:.4e} vs predicted "
        f"≈ {C1_PREDICTED:.1f} (asserted > 1e-2); source_multipliers {c1.source_multipliers!r}; "
        f"max_balance={c1.max_balance_residual:.3e}"
    )
    assert c1.source_multipliers == []
    assert e_c1 > INTERIOR_BAND, f"(c1) did not miss: {e_c1:.4e}"
    _say(
        f"WF15S W1a leg {time.perf_counter() - t00:.1f}s; summed ru_maxrss {_rss_gib(comm):.2f} GiB"
    )


# --------------------------------------------------------------- W1b, wall


def _wall_study(rho_c_wall):
    geo = ThermalBallGeometry(A_, RES, wall=ContainerWall(B_ - A_, K2, rho_c_wall))
    st = Study(geo, None, UniformHeatSource(Q1))
    return st.temperature(
        RobinCooling(K1, H, T_INF), transient=TransientSpec(T_STAR / N_CN, N_CN, RC1)
    )


def test_wf15_wall_reproduces_mat18_2b():
    comm = MPI.COMM_WORLD
    t00 = time.perf_counter()
    # `MAT-18`'s own fine rung: census, the 400 000-cell guard (asserted before
    # any solve), DG0 coefficients and the two-layer comparand.
    msh, ct, ft, nglob, k_f, rc_f, q_f, theta = _fine()
    assert nglob <= CELL_GUARD
    t0 = time.perf_counter()
    api = _wall_study(RC2)
    t_api = time.perf_counter() - t0
    direct, t_dir = _run18(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=N_CN, theta=0.5, startup=2)

    # (UT)
    ca, cb, geo = _same_mesh(api.mesh, msh, comm)
    dk = _dg0_owned_diff(api.k, k_f, comm)
    drc = _dg0_owned_diff(api.rho_c, rc_f, comm)
    dq = _dg0_owned_diff(api.q, q_f, comm)
    _say(
        f"WF15S (UT) wall: cells {ca} / {cb} (record {MAT18_T1_CELLS}); geometry equal {geo}; "
        f"DG0 k max|d| {dk[0]!r} ({dk[1]} differ), rho_c {drc[0]!r} ({drc[1]}), "
        f"q {dq[0]!r} ({dq[1]}) of {dk[2]} owned"
    )
    # (c3) the comparator, planted: one owned DG0 entry on rank 0.
    k_p = fem.Function(api.k.function_space)
    k_p.x.array[:] = api.k.x.array
    sp = float(np.spacing(np.real(k_p.x.array[0]))) if comm.rank == 0 else 0.0
    if comm.rank == 0:
        k_p.x.array[0] = k_p.x.array[0] + sp
    sp = comm.bcast(sp, root=0)
    dpl = _dg0_owned_diff(k_p, k_f, comm)
    _say(f"WF15S (c3) planted spacing {sp!r}: read back max|d| {dpl[0]!r}, differing {dpl[1]}")

    T_api = _on(msh, api.temperature)
    it = _theta_owned_rel(T_api, direct.temperature, comm)
    e_all, e_in, e_sh = _errors(msh, ct, T_api, theta)
    _, e_dir, _ = _errors(msh, ct, direct.temperature, theta)
    _say(
        f"WF15S (IT) wall: owned max|dθ|/max|θ| = {it:.3e} (band {L0_BAND:.0e}, predicted <= 1e-14)"
    )
    _say(
        f"WF15S (G) (T1) CN+Rannacher N={N_CN} res={RES} cells={ca}: interior_rel_L2={e_in:.10e} "
        f"(band 1e-2) whole={e_all:.4e} shell={e_sh:.4e} direct={e_dir:.10e}; API "
        f"max_balance={api.max_balance_residual:.3e} direct={direct.max_balance_residual:.3e}; "
        f"API solve+mesh {t_api:.1f}s, direct solve {t_dir:.1f}s"
    )
    _say(f"WF15S (R) (T1): {e_in:.4e} vs record {MAT18_T1_PRINT}")

    # (c2) the wall's ρc dropped through the API: ContainerWall(rho_c=RC1).
    c2 = _wall_study(RC1)
    _, e_c2, _ = _errors(msh, ct, _on(msh, c2.temperature), theta)
    _say(
        f"WF15S (c2) wall rho_c = RC1: interior miss measured {e_c2:.4e} vs prior same-fixture "
        f"{MAT18_C1_UNIFORM_RHOC:.4e} (asserted > 1e-2); max_balance={c2.max_balance_residual:.3e}"
    )
    _say(
        f"WF15S W1b leg {time.perf_counter() - t00:.1f}s; summed ru_maxrss {_rss_gib(comm):.2f} GiB"
    )

    assert ca == cb == MAT18_T1_CELLS, "(UT) cells"
    assert geo, "(UT) geometry"
    assert dk[0] == 0.0 and drc[0] == 0.0 and dq[0] == 0.0, "(UT) DG0 coefficients"
    assert dpl[0] == sp and dpl[1] == 1, "(c3) planted spacing"
    assert it <= L0_BAND, f"(IT) {it:.3e}"
    assert api.max_balance_residual <= BALANCE_BAND, "(G) balance"
    assert c2.max_balance_residual <= BALANCE_BAND, "(G) balance on (c2)"
    assert e_in <= INTERIOR_BAND, f"(G) (T1) {e_in:.4e}"
    assert f"{e_in:.4e}" == MAT18_T1_PRINT, f"(R) {e_in:.4e} != {MAT18_T1_PRINT}"
    assert e_c2 > INTERIOR_BAND, f"(c2) did not miss: {e_c2:.4e}"
