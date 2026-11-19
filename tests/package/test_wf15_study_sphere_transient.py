"""`WF-15` step 15 (§9 item 862, §10 U131): `MAT-12`'s SAR-driven transient through ``Study.temperature``.

On `MAT-4` step 1's sphere (EM degree 1 only — `MAT-12`'s records are degree 1;
the Study's default degree 2 is not claimed) one EM solve drives `MAT-12`'s
θ-scheme transient of the Robin-cooled gel ball through
``Study.temperature(RobinCooling, transient=TransientSpec(...))``.  T∞ = T(0) = 0
as `MAT-12` runs it (θ = T).  Brief: docs/planning/chunks/WF-15.md, "§9 item 862".

Gates (bands imported, none set):
(UT)    the Study's submesh has `MAT-10`'s 19 496 cells;
(I-dir) owned max |θ_Study − θ_direct| / max |θ_direct| ≤ ``L0_BAND`` on A
        (CN N = 40 to t*) and S (one BE step Δt₁), θ_direct = `MAT-12`'s own
        ``_cn`` / ``_be_one`` on the Study's submesh, facets and source;
(A_t)   interior relative L2 against the series ≤ ``A_BAND``; (R-A_t) its
        ``.4e`` equals `MAT-12`'s record;
(S2)    the volume-mean calorimetric identity ≤ ``S2_BAND``; (R-S2) its ``.4e``
        equals `MAT-12`'s record;
(E)     per-step balance ≤ ``E_BAND`` on A and S;
(K)     no solve: ``duty`` / ``initial_gradient_k_per_m`` / a transient probe on
        the sphere raise ``ValueError``; a port Study's transient temperature
        raises ``TypeError``; ``result`` stays None.
Control (asserted, backed by `MAT-12` control (b)): (c-qs) the steady θ standing
in for θ(t*) misses ``A_BAND``.

`MAT-12`'s records (docs/testing/logs/20260930T110230Z_MAT-12-W1.log.gz):
``:159`` sub = 19496; ``:162`` (A_t) 3.4708e-02; ``:165`` (S2) |ratio − 1| =
4.9994e-05; ``:169`` control (b) 1.0108e+00; ``:171`` EM 15.4 s, total 30.5 s.

Solve nodes run with ``FEM_EM_WF15_S15=1`` (complex build, ``-n 2``).

`WF-15` step 17 (§9 item 877, §10 U135): a ``DutyCycle`` on the EM sphere at EM
degree 1, behind ``FEM_EM_WF15_S17=1``.  Gates (Once), (UT), (χ), (S0), (Δ) ≤
``L_BAND`` (`MAT-17` step 3's Duhamel superposition over Heaviside runs), (D′)
≤ ``D_BAND`` (χ ≡ 1 reproduces the continuous run), (E) on 13 runs, (R-A_t) the
D = 1 run reproduces `MAT-12`'s 3.4708e-02; controls (Δ-ctl) and (c-sched) miss
their gate's band.  The D·P1 fraction and the closed-form duty comparison are
printed records, not gates (brief: docs/planning/chunks/WF-15.md, "§9 item 877").

`WF-15` step 16 (§9 item 863, §10 U132): the same run at the Study's **default
EM degree 2** (``_study()``, no ``degree=``), behind ``FEM_EM_WF15_S16=1``.
Gates (Deg), (Once), (UT), (I-dir), (A_t) ≤ ``A_BAND``, (S2) ≤ ``S2_BAND``, (E),
control (c-qs); no ``.4e`` equality — the degree-2 (A_t) / (S2) prints are
registered as records (brief: docs/planning/chunks/WF-15.md, "§9 item 863").

`WF-15` step 18 (§9 item 904, §10 U144): step 17's ``DutyCycle`` run at the
Study's **default EM degree 2**, behind ``FEM_EM_WF15_S18=1``.  Gates (Deg),
(Once), (UT), (χ), (S0), (Δ) ≤ ``L_BAND``, (D′) ≤ ``D_BAND``, (E) on 13 runs,
(A_t-D1) ≤ ``A_BAND`` (its ``.4e`` printed beside U132's degree-2 record, not
pinned); controls (Δ-ctl) and (c-sched) miss their gate's band (brief:
docs/planning/chunks/WF-15.md, "§9 item 904").
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import (
    DutyCycle,
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    RobinCooling,
    TransientSpec,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.probe import ProbeBall
from fem_em_solver.utils.analytical import robin_sphere_duty_cycle_temperature
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import _material_057, _study, _summed_rss_gib
from tests.package.test_wf15_study_temperature import MAT10_SUB_CELLS, _space_degree
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import (
    A_BAND,
    ROBIN_TAG,
    G,
    H,
    K,
    R,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _say,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    DT_1,
    E_BAND,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_STAR,
    _be_one,
    _cn,
    _exact_theta,
    _mean,
    _whole_rel,
)
from tests.validation.test_mat17_duty_cycle_coil_driven import D_BAND, L_BAND
from tests.validation.test_mat17_duty_cycle_fem import _census

S15_ENV = "FEM_EM_WF15_S15"
needs_s15 = pytest.mark.skipif(
    os.environ.get(S15_ENV) != "1", reason=f"opt-in solve: set {S15_ENV}=1"
)
S16_ENV = "FEM_EM_WF15_S16"
needs_s16 = pytest.mark.skipif(
    os.environ.get(S16_ENV) != "1", reason=f"opt-in solve: set {S16_ENV}=1"
)
S17_ENV = "FEM_EM_WF15_S17"
needs_s17 = pytest.mark.skipif(
    os.environ.get(S17_ENV) != "1", reason=f"opt-in solve: set {S17_ENV}=1"
)
S18_ENV = "FEM_EM_WF15_S18"
needs_s18 = pytest.mark.skipif(
    os.environ.get(S18_ENV) != "1", reason=f"opt-in solve: set {S18_ENV}=1"
)

# `MAT-12`'s printed digits (20260930T110230Z_MAT-12-W1.log.gz:162, :165).
MAT12_A_T_PRINT = "3.4708e-02"
MAT12_S2_PRINT = "4.9994e-05"


def _owned_rel(T_a, T_b, comm):
    """Owned max |T_a − T_b| / max |T_b|, reduced (P1, bs 1)."""
    n_a = T_a.function_space.dofmap.index_map.size_local
    n_b = T_b.function_space.dofmap.index_map.size_local
    same = bool(comm.allreduce(bool(n_a == n_b), op=MPI.LAND))
    a = np.asarray(T_a.x.array[:n_a], dtype=np.complex128)
    b = np.asarray(T_b.x.array[:n_b], dtype=np.complex128)
    if not same:
        return float("inf"), False
    dmax = comm.allreduce(float(np.max(np.abs(a - b), initial=0.0)), op=MPI.MAX)
    nmax = comm.allreduce(float(np.max(np.abs(b), initial=0.0)), op=MPI.MAX)
    return dmax / nmax, True


def test_wf15_s15_sphere_transient_refusals():
    """(K), no solve: duty, an initial gradient and a transient probe refused on the sphere."""
    cool = RobinCooling(K, H, 0.0, source_scale=G)
    spec = TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)
    study = _study(degree=1)
    # `WF-15` step 17 (§9 item 877): the duty pin moved — duty is gated on the
    # sphere now; without a transient it is a ValueError, a non-DutyCycle a TypeError.
    with pytest.raises(ValueError, match="duty needs a transient"):
        study.temperature(cool, duty=DutyCycle(period_s=2.0, duty=0.5))
    with pytest.raises(TypeError):
        study.temperature(cool, transient=spec, duty=0.5)
    with pytest.raises(ValueError, match="initial_gradient_k_per_m"):
        study.temperature(
            cool,
            transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C, initial_gradient_k_per_m=1.0),
        )
    with pytest.raises(ValueError, match="transient probes are not defined on the EM sphere"):
        study.probe([ProbeBall((0.0, 0.0, 0.0), 1e-3)], cool, transient=spec)
    port = _port_study()
    with pytest.raises(TypeError):
        port.temperature(cool, transient=spec)
    assert study.result is None and port.result is None, "(K) a refusal built a mesh"
    _say(
        "WF15S15 (K) duty / initial gradient / transient probe ValueError on the sphere; "
        "port Study TypeError; result None on both"
    )


def _port_study():
    from fem_em_solver import Study

    return Study(
        FourPortBirdcageGeometry.f_small(), _material_057(), LumpedPortDrive(), frequency_hz=64e6
    )


@complex_only
@needs_s15
def test_wf15_s15_sphere_transient_reproduces_mat12():
    """W1: (UT), (I-dir), (A_t) / (R-A_t), (S2) / (R-S2), (E), (c-qs); one EM solve."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    study = _study(degree=1)
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    cool = RobinCooling(K, H, 0.0, source_scale=G)

    t1 = time.perf_counter()
    thA = study.temperature(cool, transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))
    thS = study.temperature(
        cool, transient=TransientSpec(DT_1, 1, RHO_C, theta=1.0, startup_half_steps=0)
    )
    t_th = time.perf_counter() - t1
    dA = _cn(thA.mesh, thA.facet_tags, thA.q)
    dS = _be_one(thS.mesh, thS.facet_tags, thS.q)
    thQ = study.temperature(cool)
    once = study.result is res

    # (UT)
    n_sub = thA.mesh.topology.index_map(3).size_global
    # (I-dir)
    idA, layA = _owned_rel(thA.temperature, dA.temperature, comm)
    idS, layS = _owned_rel(thS.temperature, dS.temperature, comm)
    _say(f"WF15S15 (UT) sub cells {n_sub} (MAT10_SUB_CELLS {MAT10_SUB_CELLS}); once {once}")
    _say(
        f"WF15S15 (I-dir) A {idA:.3e}  S {idS:.3e}  (band {L0_BAND:.0e}); layouts "
        f"{layA} / {layS}; scheme A dt {thA.transient.dt_s:.6f} s n {thA.transient.n_steps} "
        f"theta {thA.transient.theta} startup {thA.transient.startup_half_steps}; "
        f"S dt {thS.transient.dt_s:.6e} s n {thS.transient.n_steps} theta "
        f"{thS.transient.theta} startup {thS.transient.startup_half_steps}"
    )

    # (A_t) / (R-A_t)
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub = thA.mesh
    chi = _interior_indicator(sub)
    theta = _exact_theta(sub, T_STAR, q_cf)
    err_A = _interior_rel(sub, chi, thA.temperature - theta, theta)
    _say(
        f"WF15S15 (A_t) interior rel L2 = {err_A:.10e} [{err_A:.4e}] (band {A_BAND}; "
        f"MAT-12 record {MAT12_A_T_PRINT}); q_cf {q_cf:.6f} W/m^3"
    )

    # (S2) / (R-S2)
    subS = thS.mesh
    one = fem.Constant(subS, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=subS))))
    qbar = _mean(subS, thS.q, vol)
    thbar = _mean(subS, thS.temperature, vol)
    ratio = RHO_C * thbar / (DT_1 * qbar)
    s2 = abs(ratio - 1.0)
    _say(
        f"WF15S15 (S2) rho c mean(theta^1)/(dt_1 mean q_h) = {ratio:.9f}  |ratio - 1| = "
        f"{s2:.4e} (band {S2_BAND:.0e}; MAT-12 record {MAT12_S2_PRINT}); volume {vol:.6e} "
        f"m^3  mean q_h {qbar:.6f} W/m^3"
    )

    # (E)
    eA, eS = thA.max_balance_residual, thS.max_balance_residual
    _say(f"WF15S15 (E) max balance A {eA:.3e}  S {eS:.3e}  (band {E_BAND:.0e})")

    # (c-qs) the quasi-static stand-in on the steady route's own submesh.
    subQ = thQ.submesh
    thetaQ = _exact_theta(subQ, T_STAR, q_cf)
    err_qs = _interior_rel(subQ, _interior_indicator(subQ), thQ.temperature - thetaQ, thetaQ)
    _say(
        f"WF15S15 (c-qs) steady theta for theta(t*): interior rel L2 = {err_qs:.4e} "
        f"(must miss {A_BAND}; MAT-12 control (b) 1.0108e+00)"
    )
    _say(
        f"WF15S15 elapsed: EM {t_em:.1f} s, Study transients {t_th:.1f} s, total "
        f"{time.perf_counter() - t0:.1f} s"
    )

    assert once, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(UT) sub cells {n_sub}"
    assert layA and layS, "(I-dir) owned layouts differ"
    assert idA <= L0_BAND and idS <= L0_BAND, f"(I-dir) A {idA:.3e} S {idS:.3e}"
    assert eA <= E_BAND and eS <= E_BAND, f"(E) A {eA:.3e} S {eS:.3e}"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > {A_BAND}"
    assert f"{err_A:.4e}" == MAT12_A_T_PRINT, f"(R-A_t) {err_A:.4e}"
    assert s2 <= S2_BAND, f"(S2) {s2:.4e} > {S2_BAND:.0e}"
    assert f"{s2:.4e}" == MAT12_S2_PRINT, f"(R-S2) {s2:.4e}"
    assert err_qs > A_BAND, f"(c-qs) {err_qs:.4e} did not miss {A_BAND}"


@complex_only
@needs_s16
def test_wf15_s16_sphere_transient_default_degree2():
    """W1 (`WF-15` step 16): (Deg), (Once), (UT), (I-dir), (A_t), (S2), (E), (c-qs) at EM degree 2."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    study = _study()  # the default degree
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    em_degree = _space_degree(res.fields.e_complex.function_space)
    _say(
        f"WF15S16 degree-2 EM solve + mesh {t_em:.1f} s; study.degree {study.degree}; "
        f"res.degree {res.degree}; EM space degree {em_degree}; summed ru_maxrss "
        f"{_summed_rss_gib(comm):.2f} GiB"
    )
    cool = RobinCooling(K, H, 0.0, source_scale=G)

    t1 = time.perf_counter()
    thA = study.temperature(cool, transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))
    thS = study.temperature(
        cool, transient=TransientSpec(DT_1, 1, RHO_C, theta=1.0, startup_half_steps=0)
    )
    t_th = time.perf_counter() - t1
    dA = _cn(thA.mesh, thA.facet_tags, thA.q)
    dS = _be_one(thS.mesh, thS.facet_tags, thS.q)
    thQ = study.temperature(cool)
    once = study.result is res

    # (UT)
    n_sub = thA.mesh.topology.index_map(3).size_global
    # (I-dir)
    idA, layA = _owned_rel(thA.temperature, dA.temperature, comm)
    idS, layS = _owned_rel(thS.temperature, dS.temperature, comm)
    _say(f"WF15S16 (UT) sub cells {n_sub} (MAT10_SUB_CELLS {MAT10_SUB_CELLS}); once {once}")
    _say(
        f"WF15S16 (I-dir) A {idA:.3e}  S {idS:.3e}  (band {L0_BAND:.0e}); layouts "
        f"{layA} / {layS}; scheme A dt {thA.transient.dt_s:.6f} s n {thA.transient.n_steps} "
        f"theta {thA.transient.theta} startup {thA.transient.startup_half_steps}; "
        f"S dt {thS.transient.dt_s:.6e} s n {thS.transient.n_steps} theta "
        f"{thS.transient.theta} startup {thS.transient.startup_half_steps}"
    )

    # (A_t) — registered as the degree-2 record; no equality.
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub = thA.mesh
    chi = _interior_indicator(sub)
    theta = _exact_theta(sub, T_STAR, q_cf)
    err_A = _interior_rel(sub, chi, thA.temperature - theta, theta)
    _say(
        f"WF15S16 (A_t) degree-2 record: interior rel L2 = {err_A:.10e} [{err_A:.4e}] "
        f"(band {A_BAND}; predicted [1e-3, 3.5e-2], printed only; degree-1 record "
        f"{MAT12_A_T_PRINT}); q_cf {q_cf:.6f} W/m^3"
    )

    # (S2) — registered as the degree-2 record.
    subS = thS.mesh
    one = fem.Constant(subS, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=subS))))
    qbar = _mean(subS, thS.q, vol)
    thbar = _mean(subS, thS.temperature, vol)
    ratio = RHO_C * thbar / (DT_1 * qbar)
    s2 = abs(ratio - 1.0)
    _say(
        f"WF15S16 (S2) degree-2 record: rho c mean(theta^1)/(dt_1 mean q_h) = {ratio:.9f}  "
        f"|ratio - 1| = {s2:.4e} (band {S2_BAND:.0e}; predicted approx 5e-5, printed only; "
        f"degree-1 record {MAT12_S2_PRINT}); volume {vol:.6e} m^3  mean q_h {qbar:.6f} W/m^3"
    )

    # (E)
    eA, eS = thA.max_balance_residual, thS.max_balance_residual
    _say(f"WF15S16 (E) max balance A {eA:.3e}  S {eS:.3e}  (band {E_BAND:.0e})")

    # (c-qs) the quasi-static stand-in on the steady route's own submesh.
    subQ = thQ.submesh
    thetaQ = _exact_theta(subQ, T_STAR, q_cf)
    err_qs = _interior_rel(subQ, _interior_indicator(subQ), thQ.temperature - thetaQ, thetaQ)
    _say(
        f"WF15S16 (c-qs) steady theta for theta(t*): interior rel L2 = {err_qs:.4e} "
        f"(must miss {A_BAND}; degree-1 1.0108e+00)"
    )
    _say(
        f"WF15S16 elapsed: EM {t_em:.1f} s, Study transients {t_th:.1f} s, total "
        f"{time.perf_counter() - t0:.1f} s"
    )

    assert study.degree == 2 and res.degree == 2, "(Deg) study degree"
    assert em_degree == 2, f"(Deg) EM space degree {em_degree}"
    assert once, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(UT) sub cells {n_sub}"
    assert layA and layS, "(I-dir) owned layouts differ"
    assert idA <= L0_BAND and idS <= L0_BAND, f"(I-dir) A {idA:.3e} S {idS:.3e}"
    assert eA <= E_BAND and eS <= E_BAND, f"(E) A {eA:.3e} S {eS:.3e}"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > {A_BAND}"
    assert s2 <= S2_BAND, f"(S2) {s2:.4e} > {S2_BAND:.0e}"
    assert err_qs > A_BAND, f"(c-qs) {err_qs:.4e} did not miss {A_BAND}"


def _on_space(V, T_src, comm):
    """θ_src's P1 values as a Function on ``V``. The Study rebuilds its submesh per
    call (same parent, same tag); dof layout and geometry are checked, reduced."""
    n_v = V.dofmap.index_map.size_local + V.dofmap.index_map.num_ghosts
    same = bool(comm.allreduce(bool(n_v == T_src.x.array.size), op=MPI.LAND))
    gx_v = V.mesh.geometry.x
    gx_s = T_src.function_space.mesh.geometry.x
    geo = bool(
        comm.allreduce(bool(gx_v.shape == gx_s.shape and np.array_equal(gx_v, gx_s)), op=MPI.LAND)
    )
    f = fem.Function(V)
    if same:
        f.x.array[:] = T_src.x.array
    return f, same and geo


def _exact_duty_theta(sub, t, q_cf, period, duty):
    """`MAT-17` step 1's closed form at t, interpolated in numpy into P2 (unique r, chunked)."""
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    f = fem.Function(V2)

    def ev(x):
        r_all = np.sqrt(np.sum(x**2, axis=0))
        r, inv = np.unique(r_all, return_inverse=True)
        out = np.empty(r.shape)
        for a in range(0, r.size, 20000):
            rr = r[a : a + 20000]
            out[a : a + 20000] = robin_sphere_duty_cycle_temperature(
                rr,
                np.full(rr.shape, t),
                radius=R,
                k=K,
                rho_c=RHO_C,
                q=q_cf,
                h=H,
                t_inf=0.0,
                period=period,
                duty=duty,
                periodic=False,
            )
        return out[inv.reshape(-1)]

    f.interpolate(ev)
    f.x.scatter_forward()
    return f


def _chi_string(chis):
    return "".join(str(int(round(float(np.real(c))))) for c in chis)


@complex_only
@needs_s17
def test_wf15_s17_sphere_duty_cycle():
    """W1 (`WF-15` step 17): (Once), (UT), (χ), (S0), (Δ), (D′), (E), (R-A_t); (Δ-ctl), (c-sched)."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    study = _study(degree=1)
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    cool = RobinCooling(K, H, 0.0, source_scale=G)
    dt = T_STAR / N_MAIN
    period = 8 * dt
    once = []

    t1 = time.perf_counter()
    thA = study.temperature(cool, transient=TransientSpec(dt, N_MAIN, RHO_C))
    once.append(study.result is res)
    thD = study.temperature(
        cool, transient=TransientSpec(dt, N_MAIN, RHO_C), duty=DutyCycle(period, 0.5)
    )
    once.append(study.result is res)
    thD1 = study.temperature(
        cool, transient=TransientSpec(dt, N_MAIN, RHO_C), duty=DutyCycle(period, 1.0)
    )
    once.append(study.result is res)
    t_study = time.perf_counter() - t1

    t2 = time.perf_counter()
    H_runs = []
    for m in range(10):
        H_runs.append(
            solve_transient_conduction(
                thD.mesh,
                thD.facet_tags,
                k=K,
                rho_c=RHO_C,
                q=thD.q,
                h=H,
                t_inf=0.0,
                t_init=0.0,
                robin_tag=ROBIN_TAG,
                dt=dt,
                n_steps=N_MAIN,
                theta=0.5,
                startup_half_steps=2,
                degree=1,
                source_schedule=lambda t, tm=m * 4 * dt: 1.0 if t > tm else 0.0,
            )
        )
    t_dir = time.perf_counter() - t2
    once.append(study.result is res)
    once_ok = all(once)

    sub = thD.mesh
    V1 = thD.temperature.function_space
    # (UT)
    n_sub = sub.topology.index_map(3).size_global
    _say(f"WF15S17 (UT) sub cells {n_sub} (MAT10_SUB_CELLS {MAT10_SUB_CELLS}); once {once}")

    # (χ) / (S0)
    chis = list(thD.source_multipliers)
    chi_binary = all(c in (0.0, 1.0) for c in chis)
    _say(
        f"WF15S17 (chi) D=0.5 P={period:.6f} s chi string {_chi_string(chis)} "
        f"(n {len(chis)}; binary {chi_binary})"
    )
    _say(f"WF15S17 (chi) D=1 chi string {_chi_string(thD1.source_multipliers)}")
    _census(thD, dt, 2, T_STAR / 2, "WF15S17 D=0.5")
    _census(thD1, dt, 2, T_STAR, "WF15S17 D=1")

    # (Δ) Duhamel: θ_D = Σ_m (−1)^m θ_{H_m}, summed in one P1 Function on thD.mesh.
    lay_H = True
    S = fem.Function(V1)
    S_ctl = fem.Function(V1)
    for m, h_run in enumerate(H_runs):
        arr = h_run.temperature.x.array
        ok = bool(comm.allreduce(bool(arr.size == S.x.array.size), op=MPI.LAND))
        lay_H = lay_H and ok
        if ok:
            S.x.array[:] += ((-1) ** m) * arr
            if m != 9:
                S_ctl.x.array[:] += ((-1) ** m) * arr
    err_L = _whole_rel(sub, thD.temperature - S, thD.temperature)
    err_L_ctl = _whole_rel(sub, thD.temperature - S_ctl, thD.temperature)
    _say(
        f"WF15S17 (Delta) whole rel L2 theta_D - sum(-1)^m theta_H_m = {err_L:.3e} "
        f"(band {L_BAND:.0e}; step 3 read 3.421e-14); layouts {lay_H}"
    )
    _say(
        f"WF15S17 (Delta-ctl) dropping m = 9: {err_L_ctl:.6f} (must miss {L_BAND:.0e}; "
        f"predicted 0.1-0.4; step 3 read 0.298210)"
    )

    # (D′) χ ≡ 1 reproduces the continuous run (θ_A carried onto thD's space).
    A_on, lay_A = _on_space(V1, thA.temperature, comm)
    D1_on, lay_D1 = _on_space(V1, thD1.temperature, comm)
    err_D = _whole_rel(sub, D1_on - A_on, A_on)
    _say(
        f"WF15S17 (D') whole rel L2 theta_D1 - theta_A = {err_D:.3e} (band {D_BAND:.0e}; "
        f"step 3 read 2.119e-15); layout+geometry A {lay_A} D1 {lay_D1}"
    )

    # (E)
    runs = [("A", thA), ("D", thD), ("D1", thD1)] + [(f"H{m}", h) for m, h in enumerate(H_runs)]
    es = [(lbl, float(r.max_balance_residual)) for lbl, r in runs]
    _say(
        "WF15S17 (E) max balance "
        + "  ".join(f"{lbl} {e:.3e}" for lbl, e in es)
        + f"  (band {E_BAND:.0e}; {len(es)} runs)"
    )

    # (R-A_t) / (c-sched) against MAT-12's continuous series θ_cf(t*).
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    chi = _interior_indicator(sub)
    theta_cf = _exact_theta(sub, T_STAR, q_cf)
    err_A1 = _interior_rel(sub, chi, D1_on - theta_cf, theta_cf)
    _say(
        f"WF15S17 (R-A_t) D=1 interior rel L2 = {err_A1:.10e} [{err_A1:.4e}] (band {A_BAND}; "
        f"MAT-12 / U131 record {MAT12_A_T_PRINT}); q_cf {q_cf:.6f} W/m^3"
    )
    err_cs = _interior_rel(sub, chi, thD.temperature - theta_cf, theta_cf)
    _say(
        f"WF15S17 (c-sched) theta_D vs continuous series: interior rel L2 = {err_cs:.4e} "
        f"(must miss {A_BAND}; predicted approx 0.5)"
    )

    # Records (printed only).
    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    frac = _mean(sub, thD.temperature, vol) / (0.5 * _mean(sub, A_on, vol))
    _say(
        f"WF15S17 record D*P1 fraction <theta_D>/(0.5 <theta_A>) = {frac:.6f} "
        f"(step 3 read 0.986437; printed, never asserted)"
    )
    t3 = time.perf_counter()
    theta_duty = _exact_duty_theta(sub, T_STAR, q_cf, period, 0.5)
    t_cfd = time.perf_counter() - t3
    err_cd = _interior_rel(sub, chi, thD.temperature - theta_duty, theta_duty)
    _say(
        f"WF15S17 record closed-form duty comparison: interior rel L2 = {err_cd:.4e} "
        f"(not a gate (the weekly's to promote)); P2 interpolation {t_cfd:.1f} s"
    )
    _say(
        f"WF15S17 elapsed: EM {t_em:.1f} s, Study transients (3) {t_study:.1f} s, direct "
        f"Heaviside runs (10) {t_dir:.1f} s, total {time.perf_counter() - t0:.1f} s"
    )

    assert once_ok, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(UT) sub cells {n_sub}"
    assert chi_binary, f"(chi) a midpoint chi is not 0/1: {_chi_string(chis)}"
    assert lay_H and lay_A and lay_D1, "layouts differ"
    assert err_L <= L_BAND, f"(Delta) {err_L:.3e} > {L_BAND:.0e}"
    assert err_D <= D_BAND, f"(D') {err_D:.3e} > {D_BAND:.0e}"
    assert all(e <= E_BAND for _, e in es), f"(E) {es}"
    assert err_A1 <= A_BAND, f"(R-A_t) {err_A1:.4e} > {A_BAND}"
    assert f"{err_A1:.4e}" == MAT12_A_T_PRINT, f"(R-A_t) {err_A1:.4e}"
    assert err_L_ctl > L_BAND, f"(Delta-ctl) {err_L_ctl:.3e} did not miss {L_BAND:.0e}"
    assert err_cs > A_BAND, f"(c-sched) {err_cs:.4e} did not miss {A_BAND}"
    print(f"WF15S17 VERDICT rank {comm.rank} PASS", flush=True)


# U135's degree-1 records (20261108T134930Z_WF-15-s17-W1.log) and U132's degree-2
# (A_t) record (20261107T180146Z_WF-15-s16-W1.log:136) — printed beside, never pinned.
U135_CHI_PRINT = "11111000011110000111100001111000011110000"
U132_A_T_DEG2_PRINT = "1.6447e-03"


@complex_only
@needs_s18
def test_wf15_s18_sphere_duty_cycle_default_degree2():
    """W1 (`WF-15` step 18): step 17 at EM degree 2 — (Deg), (Once), (UT), (χ), (S0), (Δ),
    (D′), (E), (A_t-D1); controls (Δ-ctl), (c-sched)."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    study = _study()  # the default degree
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    em_degree = _space_degree(res.fields.e_complex.function_space)
    _say(
        f"WF15S18 degree-2 EM solve + mesh {t_em:.1f} s; study.degree {study.degree}; "
        f"res.degree {res.degree}; EM space degree {em_degree}; summed ru_maxrss "
        f"{_summed_rss_gib(comm):.2f} GiB"
    )
    cool = RobinCooling(K, H, 0.0, source_scale=G)
    dt = T_STAR / N_MAIN
    period = 8 * dt
    once = []

    t1 = time.perf_counter()
    thA = study.temperature(cool, transient=TransientSpec(dt, N_MAIN, RHO_C))
    once.append(study.result is res)
    thD = study.temperature(
        cool, transient=TransientSpec(dt, N_MAIN, RHO_C), duty=DutyCycle(period, 0.5)
    )
    once.append(study.result is res)
    thD1 = study.temperature(
        cool, transient=TransientSpec(dt, N_MAIN, RHO_C), duty=DutyCycle(period, 1.0)
    )
    once.append(study.result is res)
    t_study = time.perf_counter() - t1

    t2 = time.perf_counter()
    H_runs = []
    for m in range(10):
        H_runs.append(
            solve_transient_conduction(
                thD.mesh,
                thD.facet_tags,
                k=K,
                rho_c=RHO_C,
                q=thD.q,
                h=H,
                t_inf=0.0,
                t_init=0.0,
                robin_tag=ROBIN_TAG,
                dt=dt,
                n_steps=N_MAIN,
                theta=0.5,
                startup_half_steps=2,
                degree=1,
                source_schedule=lambda t, tm=m * 4 * dt: 1.0 if t > tm else 0.0,
            )
        )
    t_dir = time.perf_counter() - t2
    once.append(study.result is res)
    once_ok = all(once)

    sub = thD.mesh
    V1 = thD.temperature.function_space
    # (UT)
    n_sub = sub.topology.index_map(3).size_global
    _say(f"WF15S18 (UT) sub cells {n_sub} (MAT10_SUB_CELLS {MAT10_SUB_CELLS}); once {once}")

    # (χ) / (S0)
    chis = list(thD.source_multipliers)
    chi_binary = all(c in (0.0, 1.0) for c in chis)
    _say(
        f"WF15S18 (chi) D=0.5 P={period:.6f} s chi string {_chi_string(chis)} "
        f"(n {len(chis)}; binary {chi_binary}); U135 degree 1 {U135_CHI_PRINT}"
    )
    _say(f"WF15S18 (chi) D=1 chi string {_chi_string(thD1.source_multipliers)}")
    _census(thD, dt, 2, T_STAR / 2, "WF15S18 D=0.5")
    _census(thD1, dt, 2, T_STAR, "WF15S18 D=1")

    # (Δ) Duhamel: θ_D = Σ_m (−1)^m θ_{H_m}, summed in one P1 Function on thD.mesh.
    lay_H = True
    S = fem.Function(V1)
    S_ctl = fem.Function(V1)
    for m, h_run in enumerate(H_runs):
        arr = h_run.temperature.x.array
        ok = bool(comm.allreduce(bool(arr.size == S.x.array.size), op=MPI.LAND))
        lay_H = lay_H and ok
        if ok:
            S.x.array[:] += ((-1) ** m) * arr
            if m != 9:
                S_ctl.x.array[:] += ((-1) ** m) * arr
    err_L = _whole_rel(sub, thD.temperature - S, thD.temperature)
    err_L_ctl = _whole_rel(sub, thD.temperature - S_ctl, thD.temperature)
    _say(
        f"WF15S18 (Delta) whole rel L2 theta_D - sum(-1)^m theta_H_m = {err_L:.3e} "
        f"(band {L_BAND:.0e}; U135 degree 1 read 6.608e-16); layouts {lay_H}"
    )
    _say(
        f"WF15S18 (Delta-ctl) dropping m = 9: {err_L_ctl:.6f} (must miss {L_BAND:.0e}; "
        f"predicted 0.1-0.4; U135 read 0.284785)"
    )

    # (D′) χ ≡ 1 reproduces the continuous run (θ_A carried onto thD's space).
    A_on, lay_A = _on_space(V1, thA.temperature, comm)
    D1_on, lay_D1 = _on_space(V1, thD1.temperature, comm)
    err_D = _whole_rel(sub, D1_on - A_on, A_on)
    _say(
        f"WF15S18 (D') whole rel L2 theta_D1 - theta_A = {err_D:.3e} (band {D_BAND:.0e}; "
        f"U135 read 1.604e-16); layout+geometry A {lay_A} D1 {lay_D1}"
    )

    # (E)
    runs = [("A", thA), ("D", thD), ("D1", thD1)] + [(f"H{m}", h) for m, h in enumerate(H_runs)]
    es = [(lbl, float(r.max_balance_residual)) for lbl, r in runs]
    _say(
        "WF15S18 (E) max balance "
        + "  ".join(f"{lbl} {e:.3e}" for lbl, e in es)
        + f"  (band {E_BAND:.0e}; {len(es)} runs)"
    )

    # (A_t-D1) / (c-sched) against MAT-12's continuous series θ_cf(t*).
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    chi = _interior_indicator(sub)
    theta_cf = _exact_theta(sub, T_STAR, q_cf)
    err_A1 = _interior_rel(sub, chi, D1_on - theta_cf, theta_cf)
    _say(
        f"WF15S18 (A_t-D1) D=1 interior rel L2 = {err_A1:.10e} [{err_A1:.4e}] (band {A_BAND}; "
        f"U132 degree-2 record {U132_A_T_DEG2_PRINT}, printed, not pinned); q_cf {q_cf:.6f} "
        f"W/m^3"
    )
    err_cs = _interior_rel(sub, chi, thD.temperature - theta_cf, theta_cf)
    _say(
        f"WF15S18 (c-sched) theta_D vs continuous series: interior rel L2 = {err_cs:.4e} "
        f"(must miss {A_BAND}; predicted approx 0.5; U135 read 4.9989e-01)"
    )

    # Records (printed only).
    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    frac = _mean(sub, thD.temperature, vol) / (0.5 * _mean(sub, A_on, vol))
    _say(
        f"WF15S18 record D*P1 fraction <theta_D>/(0.5 <theta_A>) = {frac:.6f} "
        f"(not a gate (the weekly's to promote); U135 degree 1 0.961466)"
    )
    t3 = time.perf_counter()
    theta_duty = _exact_duty_theta(sub, T_STAR, q_cf, period, 0.5)
    t_cfd = time.perf_counter() - t3
    err_cd = _interior_rel(sub, chi, thD.temperature - theta_duty, theta_duty)
    _say(
        f"WF15S18 record closed-form duty comparison: interior rel L2 = {err_cd:.4e} "
        f"(not a gate (the weekly's to promote); U135 degree 1 3.4400e-02); P2 "
        f"interpolation {t_cfd:.1f} s"
    )
    _say(f"WF15S18 record summed ru_maxrss {_summed_rss_gib(comm):.2f} GiB (not a gate)")
    _say(
        f"WF15S18 elapsed: EM {t_em:.1f} s, Study transients (3) {t_study:.1f} s, direct "
        f"Heaviside runs (10) {t_dir:.1f} s, total {time.perf_counter() - t0:.1f} s"
    )

    assert study.degree == 2 and res.degree == 2, "(Deg) study degree"
    assert em_degree == 2, f"(Deg) EM space degree {em_degree}"
    assert once_ok, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(UT) sub cells {n_sub}"
    assert chi_binary, f"(chi) a midpoint chi is not 0/1: {_chi_string(chis)}"
    assert lay_H and lay_A and lay_D1, "layouts differ"
    assert err_L <= L_BAND, f"(Delta) {err_L:.3e} > {L_BAND:.0e}"
    assert err_D <= D_BAND, f"(D') {err_D:.3e} > {D_BAND:.0e}"
    assert all(e <= E_BAND for _, e in es), f"(E) {es}"
    assert err_A1 <= A_BAND, f"(A_t-D1) {err_A1:.4e} > {A_BAND}"
    assert err_L_ctl > L_BAND, f"(Delta-ctl) {err_L_ctl:.3e} did not miss {L_BAND:.0e}"
    assert err_cs > A_BAND, f"(c-sched) {err_cs:.4e} did not miss {A_BAND}"
    print(f"WF15S18 VERDICT rank {comm.rank} PASS", flush=True)
