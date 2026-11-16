"""`MAT-5` step 2 — the EM-coupled σ(T) fixed point on `MAT-4`'s sphere (§9 item 105).

`MAT-10`'s mesh and submesh (R = 0.01 m, εᵣ 78, 64 MHz, h_sphere = R/10,
h_far = R/5, complex build).  Each Picard iterate solves the curl-curl problem
with a **per-cell DG0 σ = σ₀[1 + αθ̄_K]** on tag 1 (θ̄_K the DG0 interpolant of
the previous thermal iterate, carried sub → parent through both dofmaps),
turns G·½σ|E|² into a DG0 heat source on the submesh, and solves the
Robin-cooled ball (h = 3 W/m²/K, Bi = 0.05).  θ₀ = 0, N = 8 fixed, no
relaxation; the loop lives here because ``thermal/steady.py``'s Picard routine
has no per-iterate hook.  Tests only — no coupled-solver API in ``src/``.

Anchor: the lumped root θ_L = 50 K at σ* = 0.57 S/m
(``tests/materials/test_mat5_step2_lumped_root.py``), profile
``robin_sphere_steady_temperature`` at q* = 45 000 W/m³.  Errors are taken on
θ = T − T∞, never on T.

Anchors (asserted):
(L0) DG0 σ ≡ 0.57 through this EM path at `MAT-10`'s h = 100, G = 2.5e7
     reproduces `MAT-10`'s θ_h to relative L2 ≤ 1e-9 in the same run;
(P′) uniform q* at h = 3: interior relative L2 ≤ 1 %;
(A′) θ₈ against the profile at q*: interior relative L2 ≤ 11.7 % — composed:
     `MAT-4`'s 10 % (``test_lossy_sphere_sar.py:321``, by value) + step 1
     (i)'s 1 % (``P_BAND``) + 0.7 % (model 0.573 %, second order, Picard);
     never tightened after reading;
(B′) ∫_sub q_h dV = G·``dissipated_power_w`` at iterate 8, ≤ 1e-10;
(C′) the surface balance at iterate 8, ≤ 1e-9;
(D′) ∫σ_h x·x dV, parent tag 1 against sub, ≤ 1e-10;
(S)  surface-mean increments alternate in sign; |δ₈|/θ₈ ≤ 1e-3.
Negative control (predicted, asserted only as "> 11.7 %"): iterate 1 is the
frozen-σ₀ solve (closed-form miss 25.04 %).
The contraction ratios δ_{n+1}/δ_n are printed records, not asserted.
"""

from __future__ import annotations

import resource
import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.core.time_harmonic import build_material_fields, build_mu_r_field
from fem_em_solver.materials.temperature import LinearTemperatureConductivity
from fem_em_solver.post.sar import mean_sar
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature

from tests.complex_mode import complex_only
from tests.materials.test_mat5_step2_lumped_root import (
    ALPHA,
    BIOT,
    G_STEP2,
    H_STEP2,
    Q_STAR,
    SIGMA_0,
    SIGMA_STAR,
    T_REF,
    THETA_STAR,
    lumped_root,
    q_cf_unit,
)
from tests.validation.test_lossy_sphere_sar import (
    EPSILON_R_SPHERE,
    FREQUENCY_HZ,
    RHO_KG_M3,
    SIGMA_HIGH,
    SPHERE_RADIUS,
    SPHERE_TAG,
    _exact_exterior_numpy,
    _solve_lossy_sphere,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    C_BAND,
    G as G_MAT10,
    H as H_MAT10,
    K,
    P_BAND,
    QDEG,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _say,
    _sphere_submesh,
    _to_sub,
)

R = SPHERE_RADIUS
N_PICARD = 8  # fixed; no relaxation

L0_BAND = 1.0e-9
MAT4_BAND = 0.10  # test_lossy_sphere_sar.py:321, by value
A2_BAND = MAT4_BAND + P_BAND + 0.007  # 11.7 %, composed; never tightened
D_BAND = 1.0e-10
S_BAND = 1.0e-3

_MUMPS = {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}


def _real(value) -> float:
    return float(np.real(value))


def _sub_to_parent(parent_fn, sub_fn, parent_ids):
    """Overwrite the tag-1 cells of a parent DG0 field with a sub DG0 field,
    cell → dof through both dofmaps (never dof = cell)."""
    pdofs = parent_fn.function_space.dofmap.list[parent_ids, 0]
    sdofs = sub_fn.function_space.dofmap.list[np.arange(parent_ids.size), 0]
    parent_fn.x.array[pdofs] = sub_fn.x.array[sdofs]
    parent_fn.x.scatter_forward()


def _volume_mean(sub, f):
    comm = sub.comm
    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = _real(_allsum(comm, one * ufl.dx(domain=sub)))
    return _real(_allsum(comm, f * ufl.dx(domain=sub))) / vol, vol


def _em_solve(msh, cell_tags, sub, sigma_sub, parent_ids, tag):
    """One curl-curl solve with the DG0 ``sigma_sub`` on tag 1, public API only.

    Dirichlet data: `MAT-4`'s exterior closed form at the allreduced volume
    mean σ̄ of ``sigma_sub``.  Zero load, own ``LinearProblem`` and prefix.
    """
    comm = msh.comm
    sigma_bar, _ = _volume_mean(sub, sigma_sub)
    material_map = {
        SPHERE_TAG: HomogeneousMaterial(sigma=sigma_bar, epsilon_r=EPSILON_R_SPHERE)
    }
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
        cell_tags=cell_tags,
        material_map=material_map,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=_exact_exterior_numpy(sigma_bar),
    )
    solver = TimeHarmonicSolver(problem, degree=1)
    sigma_f, eps_f = build_material_fields(
        msh, problem.material, cell_tags=cell_tags, material_map=material_map
    )
    mu_f = build_mu_r_field(
        msh, problem.material, cell_tags=cell_tags, material_map=material_map
    )
    _sub_to_parent(sigma_f, sigma_sub, parent_ids)

    a = solver.bilinear_form(sigma_f, eps_f, mu_r_field=mu_f)
    bcs, _mode, _ndofs = solver.build_boundary_conditions()
    v = ufl.TestFunction(solver.function_space())
    zero = fem.Constant(msh, np.zeros(3, dtype=PETSc.ScalarType))
    lin = LinearProblem(
        a,
        ufl.inner(zero, v) * ufl.dx,
        bcs=bcs,
        petsc_options=dict(_MUMPS),
        petsc_options_prefix=f"mat5s2_em_{tag}_",
    )
    E = lin.solve()
    if isinstance(E, tuple):
        E = E[0]
    E.x.scatter_forward()
    p_diss = mean_sar(
        E, sigma=sigma_f, rho=RHO_KG_M3, cell_tags=cell_tags,
        subdomain_ids=SPHERE_TAG, comm=comm,
    )["dissipated_power_w"]
    return E, sigma_f, sigma_bar, p_diss


def _heat_source(msh, cell_tags, sub, parent_ids, E, sigma_f, gain):
    """DG0 cell average of gain·½σ|E|² on the parent, transferred to the sub."""
    Wp, q_par = _cell_average_dg0(msh, cell_tags, gain * 0.5 * sigma_f * ufl.inner(E, E))
    return _to_sub(Wp, q_par, sub, parent_ids)


def _rel_l2(sub, diff, ref):
    comm = sub.comm
    dx = ufl.dx(domain=sub, metadata=QDEG)
    e = _allsum(comm, ufl.inner(diff, diff) * dx)
    n = _allsum(comm, ufl.inner(ref, ref) * dx)
    return float(np.sqrt(_real(e) / _real(n)))


@complex_only
def test_mat5_step2_em_coupled_fixed_point():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()

    # Mesh and reference: MAT-10's EM solve, unchanged.
    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    t_ref_em = time.perf_counter() - t0

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    chi = _interior_indicator(sub)
    n_parent = msh.topology.index_map(3).size_global
    n_sub = sub.topology.index_map(3).size_global
    _say(f"\nMAT5s2 cells: parent = {n_parent}  sub = {n_sub}  (MAT-10: 74020 / 19496)")
    _say(f"MAT5s2 fixture: sigma0 = {SIGMA_0}  alpha = {ALPHA}  h = {H_STEP2}  Bi = {BIOT:.4f}  "
         f"G = {G_STEP2:.6e}  q* = {Q_STAR:.4f}  theta* = {THETA_STAR}")

    Ws = fem.functionspace(sub, ("DG", 0))
    law = LinearTemperatureConductivity(sigma_ref=SIGMA_0, alpha_per_k=ALPHA, t_ref=T_REF)
    x = ufl.SpatialCoordinate(sub)
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    area = _real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ds))

    # ---- (L0) DG0 sigma == 0.57 through this EM path against MAT-10's path.
    sig_L0 = fem.Function(Ws)
    sig_L0.x.array[:] = SIGMA_STAR
    sig_L0.x.scatter_forward()
    t1 = time.perf_counter()
    E_L0, sigf_L0, sbar_L0, _p_L0 = _em_solve(msh, cell_tags, sub, sig_L0, parent_ids, "L0")
    t_L0 = time.perf_counter() - t1
    q_L0 = _heat_source(msh, cell_tags, sub, parent_ids, E_L0, sigf_L0, G_MAT10)
    q_m10 = _heat_source(msh, cell_tags, sub, parent_ids,
                         run["e_complex"], run["sigma_field"], G_MAT10)
    res_L0 = solve_steady_conduction(sub, ft, k=K, q=q_L0, h=H_MAT10, t_inf=T_INF,
                                     robin_tag=ROBIN_TAG, degree=1)
    res_m10 = solve_steady_conduction(sub, ft, k=K, q=q_m10, h=H_MAT10, t_inf=T_INF,
                                      robin_tag=ROBIN_TAG, degree=1)
    err_L0 = _rel_l2(sub, res_L0.temperature - res_m10.temperature,
                     res_m10.temperature - T_INF)
    err_L0_int = _interior_rel(sub, chi, res_L0.temperature - res_m10.temperature,
                               res_m10.temperature - T_INF)
    _say(f"MAT5s2 (L0) sigma_bar = {sbar_L0:.15f}  theta_h[DG0 path] vs MAT-10 path: "
         f"rel L2 = {err_L0:.3e}  interior = {err_L0_int:.3e}  (bound {L0_BAND:.0e})  "
         f"EM {t_L0:.1f} s (reference solve {t_ref_em:.1f} s incl. mesh)")
    # (L0) red => the DG0-sigma path is not the gated operator: stop here.
    assert err_L0 <= L0_BAND, f"(L0) {err_L0:.3e} > {L0_BAND:.0e}"

    # ---- Profile at q*, built from the named closed form: theta = a − b r².
    th_c = float(robin_sphere_steady_temperature(
        0.0, radius=R, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF)) - T_INF
    th_R = float(robin_sphere_steady_temperature(
        R, radius=R, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF)) - T_INF
    theta_cf = th_c - (th_c - th_R) / R**2 * ufl.inner(x, x)
    _say(f"MAT5s2 profile at q*: theta(0) = {th_c:.6f} K  theta(R) = {th_R:.6f} K  "
         f"(theta_L = {THETA_STAR})")

    # ---- (P') uniform q* at h = 3.
    res_P = solve_steady_conduction(sub, ft, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_P = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta_cf, theta_cf)
    _say(f"MAT5s2 (P') uniform q*: interior rel L2 = {err_P:.4e}  (bound {P_BAND:.0e})")

    # ---- Picard, theta_0 = 0, N = 8 fixed, no relaxation.
    T_prev = None
    th_s_hist = [0.0]
    sbar_hist, em_secs, err_hist = [], [], []
    last = {}
    for n in range(1, N_PICARD + 1):
        sig_sub = fem.Function(Ws)
        if T_prev is None:
            t_cell = np.full(sig_sub.x.array.shape, T_INF, dtype=float)
        else:
            t_dg0 = fem.Function(Ws)
            t_dg0.interpolate(T_prev)
            t_dg0.x.scatter_forward()
            t_cell = np.real(t_dg0.x.array)
        sig_sub.x.array[:] = law.sigma(t_cell)
        sig_sub.x.scatter_forward()

        t1 = time.perf_counter()
        E_n, sigf_n, sbar_n, p_n = _em_solve(msh, cell_tags, sub, sig_sub, parent_ids, f"it{n}")
        em_secs.append(time.perf_counter() - t1)
        q_n = _heat_source(msh, cell_tags, sub, parent_ids, E_n, sigf_n, G_STEP2)
        res_n = solve_steady_conduction(sub, ft, k=K, q=q_n, h=H_STEP2, t_inf=T_INF,
                                        robin_tag=ROBIN_TAG, degree=1)
        th_s = _real(_allsum(comm, (res_n.temperature - T_INF) * ds)) / area
        err_n = _interior_rel(sub, chi, (res_n.temperature - T_INF) - theta_cf, theta_cf)
        th_s_hist.append(th_s)
        sbar_hist.append(sbar_n)
        err_hist.append(err_n)
        _say(f"MAT5s2 iterate {n}: sigma_bar_in = {sbar_n:.9f}  theta_s = {th_s:.9f} K  "
             f"delta = {th_s - th_s_hist[-2]:+.6e}  interior rel L2 vs profile = {err_n:.4e}  "
             f"EM {em_secs[-1]:.1f} s")
        T_prev = res_n.temperature
        last = dict(E=E_n, sigf=sigf_n, sig_sub=sig_sub, sbar=sbar_n, p=p_n, q=q_n, res=res_n)

    deltas = [th_s_hist[i] - th_s_hist[i - 1] for i in range(1, N_PICARD + 1)]
    ratios = [deltas[i] / deltas[i - 1] for i in range(1, N_PICARD)]
    _say("MAT5s2 (S) theta_s,n   = " + "  ".join(f"{v:.6f}" for v in th_s_hist[1:]))
    _say("MAT5s2 (S) sigma_bar_n = " + "  ".join(f"{v:.6f}" for v in sbar_hist))
    _say("MAT5s2 (S) delta_n     = " + "  ".join(f"{v:+.4e}" for v in deltas))
    _say("MAT5s2 (S) delta_{n+1}/delta_n = " + "  ".join(f"{v:+.5f}" for v in ratios)
         + "  (record; predicted -0.30 to -0.31)")
    alternate = all(deltas[i] * deltas[i - 1] < 0.0 for i in range(1, N_PICARD))
    s_rel = abs(deltas[-1]) / abs(th_s_hist[-1])
    _say(f"MAT5s2 (S) alternate = {alternate}  |delta_8|/theta_8 = {s_rel:.3e}  "
         f"(bound {S_BAND:.0e}, predicted 3.0e-4)")

    # ---- (A') and its legs.
    res8 = last["res"]
    err_A = err_hist[-1]
    leg_offset = _interior_rel(sub, chi, res8.temperature - res_P.temperature, theta_cf)
    _say(f"MAT5s2 (A') theta_8 interior rel L2 = {err_A:.4e}  (bound {A2_BAND:.3f}, "
         f"predicted ≈ 2.3e-2)")
    _say(f"MAT5s2 (A') legs: theta_8 - theta_h[q*] / theta = {leg_offset:.4e}; "
         f"(P') leg = {err_P:.4e}")

    # ---- (B') transfer at iterate 8.
    int_q = _real(_allsum(comm, last["q"] * ufl.dx(domain=sub)))
    p_par = G_STEP2 * last["p"]
    err_B = abs(int_q - p_par) / abs(p_par)
    _say(f"MAT5s2 (B') int_sub q_h = {int_q:.15e}  G*P_diss = {p_par:.15e}  "
         f"rel = {err_B:.3e}  (bound {B_BAND:.0e})")

    # ---- (C') balance at iterate 8.
    err_C = abs(res8.convected_power - res8.source_power) / abs(res8.source_power)
    _say(f"MAT5s2 (C') conv = {res8.convected_power:.15e}  src = {res8.source_power:.15e}  "
         f"rel = {err_C:.3e}  (bound {C_BAND:.0e})")

    # ---- (D') sigma transfer moment at iterate 8: parent tag 1 against sub.
    xp = ufl.SpatialCoordinate(msh)
    dx1 = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(SPHERE_TAG)
    m_par = _real(_allsum(comm, last["sigf"] * ufl.inner(xp, xp) * dx1))
    m_sub = _real(_allsum(comm, last["sig_sub"] * ufl.inner(x, x) * ufl.dx(domain=sub)))
    err_D = abs(m_par - m_sub) / abs(m_sub)
    _say(f"MAT5s2 (D') int sigma x.x: parent = {m_par:.15e}  sub = {m_sub:.15e}  "
         f"rel = {err_D:.3e}  (bound {D_BAND:.0e})")

    # ---- Records.
    th_s_h = th_s_hist[-1]
    root_c = lumped_root(BIOT / 5.0)
    q_bar, vol = _volume_mean(sub, last["q"])
    q_cf_bar = G_STEP2 * q_cf_unit(last["sbar"])
    # sigma of the final state (from theta_8), a record.
    t8 = fem.Function(Ws)
    t8.interpolate(res8.temperature)
    t8.x.scatter_forward()
    sig8 = fem.Function(Ws)
    sig8.x.array[:] = law.sigma(np.real(t8.x.array))
    sig8.x.scatter_forward()
    sbar8, _ = _volume_mean(sub, sig8)
    _say(f"MAT5s2 record theta_s,h = {th_s_h:.6f} K  vs theta_L {THETA_STAR}: "
         f"{th_s_h / THETA_STAR - 1.0:+.4%}  vs c = Bi/5 root {root_c:.6f}: "
         f"{th_s_h / root_c - 1.0:+.4%}")
    _say(f"MAT5s2 record q_bar_h = {q_bar:.4f} W/m^3  q_cf(sigma_bar = {last['sbar']:.6f}) = "
         f"{q_cf_bar:.4f}  q_bar_h/q_cf - 1 = {q_bar / q_cf_bar - 1.0:+.4e}  "
         f"sigma_bar[theta_8] = {sbar8:.6f}  V_sub = {vol:.6e}")
    ctrl = err_hist[0]
    _say(f"MAT5s2 control (frozen sigma0, iterate 1): interior rel L2 = {ctrl:.4e}  "
         f"(must miss {A2_BAND:.3f}; closed-form miss "
         f"{q_cf_unit(SIGMA_0) / q_cf_unit(SIGMA_STAR) - 1.0:.4%}, FEM predicted ≈ 0.29)")
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM)
    _say("MAT5s2 EM seconds per iterate: " + "  ".join(f"{s:.1f}" for s in em_secs)
         + f"  (L0 {t_L0:.1f}, reference {t_ref_em:.1f});  total "
         f"{time.perf_counter() - t0:.1f} s;  peak RSS summed over ranks "
         f"{rss / 1024.0**2:.2f} GiB")

    assert err_P <= P_BAND, f"(P') {err_P:.4e} > 1 %"
    assert err_B <= B_BAND, f"(B') {err_B:.3e} > {B_BAND:.0e}"
    assert err_C <= C_BAND, f"(C') {err_C:.3e} > {C_BAND:.0e}"
    assert err_D <= D_BAND, f"(D') {err_D:.3e} > {D_BAND:.0e}"
    assert alternate, f"(S) increments do not alternate: {deltas}"
    assert s_rel <= S_BAND, f"(S) |delta_8|/theta_8 = {s_rel:.3e} > {S_BAND:.0e}"
    assert err_A <= A2_BAND, f"(A') {err_A:.4e} > 11.7 % (composed band)"
    assert ctrl > A2_BAND, f"control (frozen sigma0) {ctrl:.4e} did not miss 11.7 %"
