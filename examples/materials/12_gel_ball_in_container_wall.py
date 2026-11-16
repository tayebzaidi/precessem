"""Example (`EX-86`): what a container wall does to the gel ball.

`mat:8` settles a bare 10 cm gel ball; `mat:4`, `mat:6`, `mat:7`, `mat:10` and
`mat:11` all solve one homogeneous medium. This one puts `MAT-18`'s 1 cm
low-conductivity container wall (a < r < b = 0.11 m, perfect contact) round the
same gel core and shows three things: the steady temperature kinks at r = a
while the heat flux stays continuous; the wall's series resistance raises the
core above the bare `MAT-9` ball; and the composite's slowest mode slows the
settling. The gated transient at t* is checked against the composite series.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The fixture,
helpers and comparands are imported from `MAT-18`'s three modules
(``test_mat18_two_layer_sphere_closed_form.py``, ``..._fem.py``,
``..._transient_fem.py``) and the bands from `MAT-13`. No gate edit.

Asserted: (S) the steady interior error and the balance at the 0.005 rung;
(T1) CN + Rannacher N = 40 at t*, interior, against the composite series;
(T2) every transient run's balance; the structural check that the gel core is
`MAT-9`'s ball; max |Im T_h| = 0 in the complex build; and two negative
controls, each asserted only as a band miss: (c-k) the wall given the gel's
conductivity, (c-rhoc) the wall given the gel's heat capacity.

Printed only, never gated: the steady profile across the interface; the
resistance budget; the flux continuity against the gradient kink; the settle
table, tau_1 and t_99 against the bare ball. Log records go in ``RECORD``.

Run it::

    ./run_examples.sh -e mat:12 -n 2 -t 420

Output (on the 0.005 mesh): ``examples/materials/paraview_output/
materials_12_gel_ball_in_container_wall_combined.xdmf`` and
``materials_12_gel_ball_in_container_wall_settle_curve.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, mesh as dmesh
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.steady import solve_steady_conduction  # noqa: E402
from fem_em_solver.utils.analytical import (  # noqa: E402
    robin_sphere_eigenvalues,
    robin_sphere_steady_temperature,
    robin_sphere_transient_temperature,
    two_layer_sphere_eigenvalues,
)

from tests.validation import test_mat9_transient_conduction_sphere as mat9  # noqa: E402
from tests.validation import test_mat18_two_layer_sphere_fem as m18a  # noqa: E402
from tests.validation import test_mat18_two_layer_sphere_transient_fem as m18b  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import (  # noqa: E402
    BALANCE_BAND,
    INTERIOR_BAND,
)
from tests.validation.test_mat18_two_layer_sphere_closed_form import (  # noqa: E402
    A_,
    B_,
    H,
    K1,
    K2,
    Q1,
    Q2,
    RC1,
    RC2,
    T_INF,
    T_R,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_12_gel_ball_in_container_wall"
RES = 0.005
N_SETTLE = 120
FO_END = 3.0
FO_PRINT = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0)
BAND_W = 2.5e-3  # interface band half-width [m]

#: Gate records, printed beside this run and never asserted.
RECORD = {
    "cells": 215542,            # 20261003T172548Z_MAT-18-step2a-W1.log:89
    "e_steady": 1.1206e-03,     # same log :89-90
    "bal_steady": 1.019e-14,    # same log :90
    "c_k": 0.1986,              # same log :110
    "e_t1": 8.2066e-04,         # 20261003T184850Z_MAT-18-step2b-W1b.log:88
    "bal_t1": 1.040e-14,        # same log :89
    "c_rc": 5.3140e-02,         # same log :94
    "theta0": 14.2550,          # 20261003T140149Z_MAT-18-step1-W2.log:89
    "theta0_bare": 12.2222,     # same log :106
    "tau1": 2.4549e+04,         # same log :106
    "tau1_bare": 19155.0,       # 20261002T051536Z_EX-74.log:39
    "t99_bare": 88118.0,        # same log :70
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _theta_of(T_h, name):
    f = fem.Function(T_h.function_space, name=name)
    f.x.array[:] = T_h.x.array - T_INF
    f.x.scatter_forward()
    return f


def _p1_of(u, V1, name):
    f = fem.Function(V1, name=name)
    f.interpolate(u)
    f.x.scatter_forward()
    return f


def _at(T, pts):
    vals, ok = evaluate_vector_field_parallel(T, pts)
    out = []
    for i in range(pts.shape[0]):
        out.append(float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else float("nan"))
    return np.array(out)


def _max_im(comm, fs):
    m = max(float(np.max(np.abs(np.imag(f.x.array)))) if f.x.array.size else 0.0 for f in fs)
    return comm.allreduce(m, op=MPI.MAX)


def _first_below(ts, deficits, level):
    for i in range(1, len(ts)):
        if deficits[i] <= level < deficits[i - 1]:
            a, b = np.log(deficits[i - 1]), np.log(deficits[i])
            return float(ts[i - 1] + (np.log(level) - a) / (b - a) * (ts[i] - ts[i - 1]))
    return float("nan")


_GLX, _GLW = np.polynomial.legendre.leggauss(64)
_RG = (0.5 * A_ * (_GLX + 1.0))[:, None]  # (64, 1) radii on [0, a]
_WG = (0.5 * A_ * _GLW)[:, None]


def _gl_mean(theta_fn):
    """3/a^3 int_0^a theta r^2 dr by 64-point Gauss-Legendre; theta_fn(r)->(64, ...)"""
    th = theta_fn(_RG)
    return 3.0 / A_**3 * np.sum(_WG * th * _RG**2, axis=0)


def _core_series(t):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = []
    for s in range(0, t.size, 50):
        out.append(_gl_mean(lambda r: m18b._two_layer(r, t[s:s + 50]) - T_INF))
    return np.concatenate(out)


def _bare_series(t):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = []
    for s in range(0, t.size, 50):
        out.append(_gl_mean(lambda r: robin_sphere_transient_temperature(
            r, t[s:s + 50], radius=A_, k=K1, rho_c=RC1, q=Q1, h=H, t_inf=T_INF,
            t_init=T_INF) - T_INF))
    return np.concatenate(out)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # Structural: the gel core is MAT-9's ball, so the wall is the only difference.
    assert A_ == mat9.R and K1 == mat9.K and Q1 == mat9.Q and RC1 == mat9.RHO_C
    assert H == mat9.H and T_INF == mat9.T_INF and T_R == mat9.T_R
    assert m18b.T_STAR == mat9.T_END

    _say(comm, "=" * 78)
    _say(comm, "EX-86 -- what a container wall does to the gel ball (MAT-18 fixture)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: interior {INTERIOR_BAND:.0e}, balance {BALANCE_BAND:.0e}; "
               f"core a = {A_} m (k {K1}, rho c {RC1:.3g}, q {Q1}), wall b = {B_} m "
               f"(k {K2}, rho c {RC2:.3g}, q {Q2}); h = {H}, T_inf = {T_INF}; "
               f"t* = {m18b.T_STAR:.1f} s; critical radius 2 k2/h = {2 * K2 / H:.4f} m < a")

    msh, ct, ft, nglob, k_f, rc_f, q_f, theta_tc = m18b._fine()
    t_build = m18a._ball(RES)[3]
    vc = m18a._vcore(msh, ct)
    vol = m18a._allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * ufl.dx(domain=msh))
    _say(comm, f"[mesh] res {RES}: {nglob} cells (record {RECORD['cells']}), build {t_build:.1f} s; "
               f"Vcore_h/Vcore-1 = {vc / m18b.V_CORE - 1:.4e}, "
               f"Vshell_h/Vshell-1 = {(vol - vc) / m18b.V_SHELL - 1:.4e}")

    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"gel core (k {K1} W/m/K, rho c {RC1:.3g})",
            2: f"container wall (k {K2} W/m/K, rho c {RC2:.3g})",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:12 -- gel ball in a 1 cm container wall (MAT-18 fixture), Robin outside",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- (S) steady, 0.005 rung --------------------------------------------------------
    sol = solve_steady_conduction(msh, ft, k=k_f, q=q_f, h=H, t_inf=T_INF,
                                  robin_tag=m18a.ROBIN_TAG, degree=1)
    theta_cf = m18a._comparand(msh, m18a._two_layer)
    e_all, e_in, e_sh = m18a._errors(msh, ct, sol.temperature, theta_cf)
    bal_s = abs(sol.convected_power - Q1 * vc) / (Q1 * vc)
    _say(comm, f"\n[steady] (S) interior rel L2 {e_in:.4e} (record {RECORD['e_steady']:.4e}); "
               f"whole {e_all:.4e}, shell {e_sh:.4e}  ASSERTED interior <= {INTERIOR_BAND:.0e}")
    _say(comm, f"[steady] (S) balance {bal_s:.3e} (record {RECORD['bal_steady']:.3e})  "
               f"ASSERTED <= {BALANCE_BAND:.0e}")

    # ---- (c-k) negative control: the wall given the gel's conductivity ----------------
    sol_k = solve_steady_conduction(msh, ft, k=K1, q=q_f, h=H, t_inf=T_INF,
                                    robin_tag=m18a.ROBIN_TAG, degree=1)
    _, ek_in, _ = m18a._errors(msh, ct, sol_k.temperature, theta_cf)
    _say(comm, f"[control] (c-k) uniform k = K1: steady interior miss {ek_in:.4f} "
               f"(record {RECORD['c_k']:.4f}); ASSERTED > {INTERIOR_BAND:.0e}")

    # ---- the profile across the interface -----------------------------------------------
    fr = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 0.97, 1.0])
    rs = np.concatenate([fr * A_, A_ + np.array([0.25, 0.5, 0.75]) * (B_ - A_), [0.999 * B_]])
    pts = np.zeros((rs.size, 3))
    pts[:, 0] = rs
    th_h = _at(sol.temperature, pts)
    th_cf = m18a._two_layer(rs) - T_INF
    th_bare = np.where(rs <= A_ * (1 + 1e-12),
                       robin_sphere_steady_temperature(np.minimum(rs, A_), radius=A_, k=K1,
                                                       q=Q1, h=H, t_inf=T_INF) - T_INF, np.nan)
    _say(comm, "\n[profile] steady theta along +x [K]:   r/a     theta_h   two-layer   rel diff   bare ball")
    for r, a_, b_, c_ in zip(rs, th_h, th_cf, th_bare):
        bare = f"{c_:10.4f}" if np.isfinite(c_) else "         -"
        _say(comm, f"[profile]   {r / A_:7.4f}   {a_:9.4f}   {b_:9.4f}   {a_ / b_ - 1:+.3e}   {bare}"
                   f"{'  (valid=False)' if not np.isfinite(a_) else ''}")
    _say(comm, f"[profile] theta_h(0) = {th_h[0]:.4f} K beside record {RECORD['theta0']} K "
               f"(predicted within 1e-3 relative)")

    # ---- the resistance budget ---------------------------------------------------------
    cf3 = m18a._two_layer(np.array([0.0, A_, B_])) - T_INF
    bare0 = float(robin_sphere_steady_temperature(0.0, radius=A_, k=K1, q=Q1, h=H, t_inf=T_INF)) - T_INF
    bare_a = float(robin_sphere_steady_temperature(A_, radius=A_, k=K1, q=Q1, h=H, t_inf=T_INF)) - T_INF
    _say(comm, "\n[budget] drops [K]        core(0->a)   wall(a->b)   film(b->inf)   centre theta(0)")
    _say(comm, f"[budget] two-layer       {cf3[0] - cf3[1]:10.3f}   {cf3[1] - cf3[2]:10.3f}   "
               f"{cf3[2]:10.3f}     {cf3[0]:.4f} (record {RECORD['theta0']})")
    _say(comm, f"[budget] FEM (0.999 b)   {th_h[0] - th_h[6]:10.3f}   {th_h[6] - th_h[-1]:10.3f}   "
               f"{th_h[-1]:10.3f}     {th_h[0]:.4f}")
    _say(comm, f"[budget] bare MAT-9 ball {bare0 - bare_a:10.3f}   {0.0:10.3f}   {bare_a:10.3f}     "
               f"{bare0:.4f} (record {RECORD['theta0_bare']})")
    _say(comm, f"[budget] the wall raises the centre by {cf3[0] - bare0:.3f} K "
               f"({(cf3[0] / bare0 - 1) * 100:.2f} %); critical radius of insulation 2 k2/h = "
               f"{2 * K2 / H:.4f} m < a = {A_} m, so any thickness of this wall raises the core")

    # ---- flux continuity against the kink ----------------------------------------------
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    nloc = imap.size_local
    ncell = nloc + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(ncell, dtype=np.int32))
    rm = np.linalg.norm(mids, axis=1)
    tags = np.zeros(ncell, dtype=int)
    tags[ct.find(1)] = 1
    tags[ct.find(2)] = 2
    Wv = fem.functionspace(msh, ("DG", 0, (3,)))
    gradf = fem.Function(Wv, name="grad_theta_h_steady")
    gradf.interpolate(fem.Expression(ufl.grad(sol.temperature), Wv.element.interpolation_points))
    gradf.x.scatter_forward()
    g = np.real(gradf.x.array).reshape(-1, 3)[:ncell]
    rhat = mids / np.where(rm > 0, rm, 1.0)[:, None]
    dr = np.sum(g * rhat, axis=1)
    kcell = np.where(tags == 1, K1, np.where(tags == 2, K2, np.nan))
    qr = -kcell * dr
    Wd = fem.functionspace(msh, ("DG", 0))
    cvol = np.real(fem.assemble_vector(fem.form(ufl.conj(ufl.TestFunction(Wd)) * ufl.dx(domain=msh))).array)[:nloc]
    cvol = np.concatenate([cvol, np.zeros(ncell - nloc)])
    band_c = (tags == 1) & (rm >= A_ - BAND_W) & (rm < A_)
    band_s = (tags == 2) & (rm > A_) & (rm <= A_ + BAND_W)
    own = np.arange(ncell) < nloc

    def _wmean(x, mask):
        m = mask & own
        s = comm.allreduce(float(np.sum(cvol[m] * x[m])), op=MPI.SUM)
        v = comm.allreduce(float(np.sum(cvol[m])), op=MPI.SUM)
        return s / v

    gc, gs = _wmean(np.abs(dr), band_c), _wmean(np.abs(dr), band_s)
    qc, qs = _wmean(qr, band_c), _wmean(qr, band_s)
    gauss = Q1 * vc / (4.0 * np.pi * A_**2)
    _say(comm, f"\n[flux] |d theta/dr| core band {gc:.4f}, wall band {gs:.4f} K/m; "
               f"shell/core {gs / gc:.3f} (predicted [2.7, 3.6]; k1/k2 = {K1 / K2:.3f})")
    _say(comm, f"[flux] q_r core band {qc:.4f}, wall band {qs:.4f} W/m^2; shell/core "
               f"{qs / qc:.3f} (predicted [0.9, 1.1]); Gauss Q1 Vcore_h/(4 pi a^2) = "
               f"{gauss:.4f} W/m^2")

    # ---- (T1) transient at t*, (c-rhoc), settle run --------------------------------------
    dt = m18b.T_STAR / m18b.N_CN
    s1, t_s1 = m18b._run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=m18b.N_CN, theta=0.5, startup=2)
    _, e1_in, _ = m18b._errors(msh, ct, s1.temperature, theta_tc)
    _say(comm, f"\n[transient] (T1) CN + Rannacher N = {m18b.N_CN} at t*: interior rel L2 "
               f"{e1_in:.4e} (record {RECORD['e_t1']:.4e})  ASSERTED <= {INTERIOR_BAND:.0e}; "
               f"balance {s1.max_balance_residual:.3e} (record {RECORD['bal_t1']:.3e}); "
               f"solve {t_s1:.1f} s")
    s_rc, t_rc = m18b._run(msh, ft, k=k_f, rho_c=m18b.RC1, q=q_f, n=m18b.N_CN, theta=0.5, startup=2)
    _, erc_in, _ = m18b._errors(msh, ct, s_rc.temperature, theta_tc)
    _say(comm, f"[control] (c-rhoc) uniform rho c = RC1: interior miss {erc_in:.4e} "
               f"(record {RECORD['c_rc']:.4e}); ASSERTED > {INTERIOR_BAND:.0e}; "
               f"balance {s_rc.max_balance_residual:.3e}")

    dt_set = FO_END * T_R / N_SETTLE
    steps_all = tuple(range(0, N_SETTLE + 1))
    sres, t_set = m18b._run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=N_SETTLE, theta=0.5, startup=2,
                            dt=dt_set, snapshot_steps=steps_all)
    balances = {"(T1) CN N=40": s1.max_balance_residual,
                "(c-rhoc) CN N=40": s_rc.max_balance_residual,
                "settle CN N=120": sres.max_balance_residual}
    worst = max(balances, key=balances.get)
    _say(comm, f"[transient] (T2) worst per-step balance {balances[worst]:.3e} on '{worst}' of "
               f"{len(balances)}  ASSERTED <= {BALANCE_BAND:.0e}; settle solve {t_set:.1f} s")

    # ---- settle table ----------------------------------------------------------------------
    chi_core = m18a._dg0(msh, ct, 1.0, 0.0)
    fem_mean = np.array([
        m18a._allsum(comm, chi_core * (sres.snapshots[i] - T_INF) * ufl.dx(domain=msh)) / vc
        for i in steps_all])
    steady_fem = m18a._allsum(comm, chi_core * (sol.temperature - T_INF) * ufl.dx(domain=msh)) / vc
    t_all = np.array(steps_all) * dt_set
    fo_all = t_all / T_R
    ser_ss = float(np.ravel(_gl_mean(lambda r: m18a._two_layer(r) - T_INF))[0])
    bare_ss = float(np.ravel(_gl_mean(lambda r: robin_sphere_steady_temperature(
        r, radius=A_, k=K1, q=Q1, h=H, t_inf=T_INF) - T_INF))[0])
    ser_mean = np.concatenate([[0.0], _core_series(t_all[1:])])
    bare_mean = np.concatenate([[0.0], _bare_series(t_all[1:])])
    _say(comm, f"\n[settle] steady core-mean theta: FEM {steady_fem:.4f} K, two-layer series "
               f"{ser_ss:.4f} K; bare ball {bare_ss:.4f} K (+{(ser_ss / bare_ss - 1) * 100:.1f} %) "
               f"(predicted 10.92 vs 8.889); CN N = {N_SETTLE}, dt = {dt_set:.1f} s")
    _say(comm, "[settle]    Fo    FEM core mean   series core mean   FEM-series   bare-ball mean   "
               "theta_h(0)   series(0)")
    for fo in FO_PRINT:
        i = int(round(fo * N_SETTLE / FO_END))
        assert abs(fo_all[i] - fo) < 1e-9
        p0 = _at(sres.snapshots[i], np.zeros((1, 3)))[0]
        s0 = float(m18b._two_layer(np.array([0.0]), np.array([t_all[i]]))[0] - T_INF)
        _say(comm, f"[settle] {fo:5.2f}  {fem_mean[i]:13.6e}   {ser_mean[i]:13.6e}   "
                   f"{fem_mean[i] - ser_mean[i]:+.3e}   {bare_mean[i]:13.6e}   {p0:9.4f}   {s0:9.4f}")

    s1_rate = float(two_layer_sphere_eigenvalues(a=A_, b=B_, k1=K1, k2=K2, rho_c1=RC1,
                                                 rho_c2=RC2, h=H, n=1)[0])
    tau1 = 1.0 / s1_rate
    lam1 = float(robin_sphere_eigenvalues(H * A_ / K1, 1)[0])
    tau1_bare = T_R / lam1**2
    _say(comm, f"[tau] composite tau_1 = 1/s_1 = {tau1:.1f} s (record {RECORD['tau1']:.4e}); bare "
               f"ball tau_1 = t_R/lambda_1^2 = {tau1_bare:.1f} s (record {RECORD['tau1_bare']:.1f}); "
               f"ratio {tau1 / tau1_bare:.3f} (predicted 1.282)")

    # t_99 on a Fo grid of step 0.01 over [1, 3] for the series; every step for the FEM
    fo_grid = np.round(np.arange(1.0, 3.0 + 1e-9, 0.01), 10)
    t_grid = fo_grid * T_R
    ser_grid = _core_series(t_grid)
    bare_grid = _bare_series(t_grid)
    t99_ser = _first_below(t_grid, (ser_ss - ser_grid) / ser_ss, 0.01)
    t99_bare = _first_below(t_grid, (bare_ss - bare_grid) / bare_ss, 0.01)
    t99_fem = _first_below(t_all, (steady_fem - fem_mean) / steady_fem, 0.01)
    _say(comm, f"[t99] core mean, 1 % deficit: two-layer series {t99_ser:.0f} s "
               f"(predicted [100000, 130000]); FEM {t99_fem:.0f} s "
               f"(rel diff to series {t99_fem / t99_ser - 1:+.2e}); bare ball series "
               f"{t99_bare:.0f} s (record {RECORD['t99_bare']:.0f}, rel {t99_bare / RECORD['t99_bare'] - 1:+.2e}); "
               f"ratio composite/bare {t99_ser / t99_bare:.3f} (predicted [1.15, 1.45])")
    d_end = m18b._rel_diff(msh, sres.temperature, sol.temperature)
    _say(comm, f"[settle] settle end (Fo 3) vs steady solve, whole-ball rel diff {d_end:.3e} "
               f"(predicted ~1e-4)")

    # ---- complex build: no imaginary part ----------------------------------------------------
    im = _max_im(comm, [sol.temperature, sol_k.temperature, s1.temperature, s_rc.temperature,
                        sres.temperature])
    _say(comm, f"[complex] max |Im T_h| over every solve = {im:.3e}  ASSERTED == 0")

    # ---- assertions --------------------------------------------------------------------------
    assert e_in <= INTERIOR_BAND, f"(S) interior {e_in:.4e}"
    assert bal_s <= BALANCE_BAND, f"(S) balance {bal_s:.3e}"
    assert e1_in <= INTERIOR_BAND, f"(T1) interior {e1_in:.4e}"
    assert balances[worst] <= BALANCE_BAND, f"(T2) {worst} {balances[worst]:.3e}"
    assert ek_in > INTERIOR_BAND, f"(c-k) did not fire: {ek_in:.4e}"
    assert erc_in > INTERIOR_BAND, f"(c-rhoc) did not fire: {erc_in:.4e}"
    assert im == 0.0, f"max |Im T_h| = {im:.3e}"

    # ---- outputs -----------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = sol.temperature.function_space
    th_steady = _theta_of(sol.temperature, "theta_h_steady")
    th_cf_p1 = _p1_of(theta_cf, V1, "theta_cf_steady")
    th_err = fem.Function(V1, name="theta_err_steady")
    th_err.x.array[:] = th_steady.x.array - th_cf_p1.x.array
    th_err.x.scatter_forward()
    th_ts = _theta_of(s1.temperature, "theta_h_tstar")
    th_ser = _p1_of(theta_tc, V1, "theta_series_tstar")
    th_terr = fem.Function(V1, name="theta_err_tstar")
    th_terr.x.array[:] = th_ts.x.array - th_ser.x.array
    th_terr.x.scatter_forward()

    def _dg0_named(arr, name):
        f = fem.Function(Wd, name=name)
        f.x.array[:ncell] = arr[:ncell]
        f.x.scatter_forward()
        return f

    band_tag = np.where(band_c, 1.0, np.where(band_s, 2.0, 0.0))
    fields = {
        "theta_h_steady": th_steady,
        "theta_cf_steady": th_cf_p1,
        "theta_err_steady": th_err,
        "theta_h_uniform_k_steady": _theta_of(sol_k.temperature, "theta_h_uniform_k_steady"),
        "theta_h_tstar": th_ts,
        "theta_series_tstar": th_ser,
        "theta_err_tstar": th_terr,
        "k": _dg0_named(np.where(tags == 1, K1, K2).astype(float), "k"),
        "rho_c": _dg0_named(np.where(tags == 1, RC1, RC2).astype(float), "rho_c"),
        "radial_gradient_steady": _dg0_named(dr, "radial_gradient_steady"),
        "radial_flux_steady": _dg0_named(qr, "radial_flux_steady"),
        "interface_band": _dg0_named(band_tag, "interface_band"),
        "interior_mask": m18a._indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * B_**2),
    }
    fields["interior_mask"].name = "interior_mask"
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, msh, ct, fields, comm=comm)
    table = OUTPUT_DIR / f"{BASENAME}_settle_curve.csv"
    if comm.rank == 0:
        with table.open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["step", "Fo", "t_s", "fem_core_mean_K", "series_core_mean_K",
                         "bare_ball_mean_K"])
            for i in steps_all:
                wr.writerow([i, f"{fo_all[i]:.6f}", f"{t_all[i]:.3f}", f"{fem_mean[i]:.9e}",
                             f"{ser_mean[i]:.9e}", f"{bare_mean[i]:.9e}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P1 and DG0 fields on the 0.005 mesh")
    _say(comm, f"[output] {table}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
