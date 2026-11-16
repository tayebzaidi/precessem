"""`MAT-21` step 3 — step 1's steady + CN anchors at Lagrange degree 2 against a
P3 comparand on m ∈ {2, 3, 4} (§9 item 443, §10 chain P28, daily 2026-10-14 —
provisional; brief `docs/planning/chunks/MAT-21.md`).

A sibling of step 1's ``test_mat21_box_fem.py`` (which has no degree knob and
is imported by step 2's gate and EX-112): this module imports step 1's helpers
and owns degree-2 copies.  It sets no band and moves none — every band is
imported by name.  Disclosed departures from step 1's helpers:

* the comparand is P3 (``CMP_DEGREE = DEGREE + 1``), not P2, so a P2 solution
  is measured against ‖u − u_h‖ rather than ‖I₂u − u_h‖;
* every P2 reading goes through the local ``_rel_p`` at quadrature degree 6
  (the imported ``_rel`` pins degree 4, exact for (P1 − P2)² only).  The
  imported ``_rel`` is used solely for the (c-deg) cross-print against step 1's
  P2-comparand number 1.3552e-2.

Asserted: (P0) step-0 constants = FEM's; (Pdeg) global dofs = (26m+1)(16m+1)(4m+1);
(S1) LSQ slope ≥ ``RATE_BAR`` (printed beside ≈ 3), m = 4 interior ≤
``INTERIOR_BAND``; (S2) Gauss ≤ ``S2_GAUSS_BAND`` per rung; (T) CN N = 40 at t*
interior ≤ ``INTERIOR_BAND``; (S-slope) ≤ ``S2_BAND``; (E) ≤ ``E_BAND``; (∞) ≤
``INF_BAND``.  Controls (separations only): (c-deg) P1 interior > P2 interior
on m = 2; (c1) 2h comparand and (c2) Dirichlet z-roots miss ``INTERIOR_BAND``.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _robin_tags, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    E_BAND,
    INF_BAND,
    RHO_C,
    S2_BAND,
)
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import (
    RATE_BAR,
    S2_GAUSS_BAND,
    _max_imag,
)
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    H,
    K,
    Q,
    T_C,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_fem import (
    BASE,
    DT_1,
    INTERP_CHUNK,
    T_STAR,
    V_BOX,
    _interior_chi,
    _rel,
    _steady_ex,
    _trans_ex,
)

DEGREE = 2
CMP_DEGREE = DEGREE + 1  # 3
QDEG_P = {"quadrature_degree": 2 * CMP_DEGREE}  # 6
RUNGS_P = (2, 3, 4)
M_T = RUNGS_P[-1]
M_CTL = RUNGS_P[0]
COMPARAND_FLOOR = 4.6708e-6  # step 0's 1 − ΠS_d at SERIES_TERMS (chunk file line 152)


def _dofs_pred(m, deg=DEGREE):
    n = [b * m * deg + 1 for b in BASE]
    return int(np.prod(n))


def _exact_p(msh, fn, deg=CMP_DEGREE):
    V = fem.functionspace(msh, ("Lagrange", deg))
    g = fem.Function(V)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = np.real(
                fn(x[0, s:s + INTERP_CHUNK], x[1, s:s + INTERP_CHUNK], x[2, s:s + INTERP_CHUNK]))
        return out

    g.interpolate(ev)
    npts = msh.comm.allreduce(int(V.dofmap.index_map.size_local), op=MPI.SUM)
    return g, npts


def _rel_p(msh, th_h, th_ex, chi=None):
    """sqrt(∫χ|θ_h − θ_ex|² / ∫χ|θ_ex|²) at quadrature degree 6, complex-safe."""
    dx = ufl.dx(domain=msh, metadata=QDEG_P)
    w = 1.0 if chi is None else chi
    d = th_h - th_ex
    e = _allsum(msh.comm, w * ufl.inner(d, d) * dx)
    n = _allsum(msh.comm, w * ufl.inner(th_ex, th_ex) * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


class _Clock3:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-3."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-3 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _rung_p(m, clk, degree):
    comm = MPI.COMM_WORLD
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(comm, [np.array([-A_BOX, -B_BOX, -C_BOX]),
                                  np.array([A_BOX, B_BOX, C_BOX])],
                           ns, cell_type=dmesh.CellType.tetrahedron)
    ft = _robin_tags(msh)
    n = msh.topology.index_map(3).size_global
    one = fem.Constant(msh, default_scalar_type(1.0))
    v_h = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    clk.lap(f"build m={m}")
    st = solve_steady_conduction(msh, ft, k=K_GEL, q=Q, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=degree)
    clk.lap(f"steady solve m={m} degree={degree}")
    return msh, ft, ns, n, v_h, st


@pytest.fixture(scope="module")
def steady_p():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock3("steady")
    out = {"rungs": [], "imag": 0.0}
    _say(f"\nMAT21-3 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; "
         f"rungs {RUNGS_P}; comparand floor 1 - prod S_d = {COMPARAND_FLOOR:.4e} (step 0, printed)")
    for m in RUNGS_P:
        msh, ft, ns, n, v_h, st = _rung_p(m, clk, DEGREE)
        ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
        th_h = st.temperature - T_INF
        conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
        gauss = abs(conv - Q * v_h) / (Q * v_h)
        ex, npts = _exact_p(msh, _steady_ex())
        clk.lap(f"steady series interpolation m={m} ({npts} P3 points)")
        e_whole = _rel_p(msh, th_h, ex)
        h_eff = (v_h / n) ** (1.0 / 3.0)
        out["imag"] = max(out["imag"], _max_imag(st.temperature))
        ndof = int(st.temperature.function_space.dofmap.index_map.size_global)
        pts_pred = _dofs_pred(m, CMP_DEGREE)
        rec = {"m": m, "n": n, "v_h": v_h, "h_eff": h_eff, "e": e_whole, "gauss": gauss,
               "ndof": ndof, "ndof_pred": _dofs_pred(m)}
        _say(f"MAT21-3 rung m = {m} {tuple(ns)}: cells {n} (6*prod N = {6 * int(np.prod(ns))}); "
             f"P2 dofs {ndof} (predicted {_dofs_pred(m)}); P3 points {npts} (predicted {pts_pred}); "
             f"V_h = {v_h:.15e} vs 8ABC {V_BOX:.15e} ({v_h / V_BOX - 1.0:+.3e}); "
             f"h_eff = {h_eff:.6e} m; whole L2 rel = {e_whole:.6e}; (S2) Gauss = {gauss:.3e} "
             f"(bound {S2_GAUSS_BAND:.0e}, predicted <= 5e-12)")
        out["rungs"].append(rec)

        if m == M_CTL:
            # (c-deg): the same rung at degree 1, the same P3 comparand.
            chi = _interior_chi(msh)
            e2_int = _rel_p(msh, th_h, ex, chi)
            msh1 = msh
            st1 = solve_steady_conduction(msh1, ft, k=K_GEL, q=Q, h=H_GEL, t_inf=T_INF,
                                          robin_tag=ROBIN_TAG, degree=1)
            clk.lap(f"(c-deg) steady solve m={m} degree=1")
            th1 = st1.temperature - T_INF
            e1_whole = _rel_p(msh, th1, ex)
            e1_int = _rel_p(msh, th1, ex, chi)
            ex_p2, _ = _exact_p(msh, _steady_ex(), deg=2)
            e1_cross = _rel(msh, th1, ex_p2)
            ndof1 = int(st1.temperature.function_space.dofmap.index_map.size_global)
            out["cdeg"] = {"e1_int": e1_int, "e2_int": e2_int, "e1_whole": e1_whole,
                           "e2_whole": e_whole}
            _say(f"MAT21-3 control (c-deg) m = {m}: P1 dofs {ndof1} (predicted {_dofs_pred(m, 1)}); "
                 f"interior P1 {e1_int:.6e} vs P2 {e2_int:.6e}, ratio e1/e2 = {e1_int / e2_int:.4f}; "
                 f"whole P1 {e1_whole:.6e} vs P2 {e_whole:.6e}, ratio {e1_whole / e_whole:.4f} "
                 f"(predicted 3-30, central 8; asserted interior ratio > 1); P1 whole vs P3 "
                 f"predicted 1.30e-2 to 1.42e-2; cross-print P1 via imported _rel vs P2 comparand "
                 f"= {e1_cross:.6e} (step 1 W1 :94 read 1.3552e-2)")
            clk.lap("control (c-deg)")

        if m != M_T:
            continue

        chi = _interior_chi(msh)
        out["interior"] = _rel_p(msh, th_h, ex, chi)
        try:
            vals, ok = evaluate_vector_field_parallel(st.temperature, np.zeros((1, 3)))
            th_c = float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
        except Exception as exc:  # printed quantity only
            _say(f"MAT21-3 centre evaluation failed: {exc!r}")
            th_c = float("nan")
        th_c_ex = float(np.real(_steady_ex()(np.zeros(1), np.zeros(1), np.zeros(1))[0]))
        _say(f"MAT21-3 (S1) finest (m = {m}) interior rel L2 = {out['interior']:.6e} "
             f"(bound {INTERIOR_BAND:.0e}, predicted 3e-5 to 5e-4); theta_h(0) = {th_c:.6f} K vs "
             f"series {th_c_ex:.6f} K (predicted 12.003947; rel {th_c / th_c_ex - 1.0:+.3e}, "
             f"predicted |rel| <= 3e-5); imported-_rel (qdeg 4) whole at m = {m}: "
             f"{_rel(msh, th_h, ex):.6e} (printed)")
        clk.lap("finest-rung interior + centre")
        ex_c1, _ = _exact_p(msh, _steady_ex(h=2.0 * H))
        out["c1"] = _rel_p(msh, th_h, ex_c1, chi)
        _say(f"MAT21-3 control (c1) comparand at 2h: steady interior miss = {out['c1']:.4e} "
             f"(predicted 0.597, range 0.59-0.60; asserted > {INTERIOR_BAND:.0e})")
        clk.lap("control (c1)")

    hs = np.array([r["h_eff"] for r in out["rungs"]])
    es = np.array([r["e"] for r in out["rungs"]])
    out["slope"] = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
    out["pairwise"] = [float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1]))
                       for i in range(len(hs) - 1)]
    _say("MAT21-3 (S1) whole L2 per rung " + ", ".join(f"{e:.4e}" for e in es)
         + " (predicted 4e-4 to 4e-3 at m = 2, then /3.4, /2.4)")
    _say(f"MAT21-3 (S1) LSQ slope log e vs log h_eff = {out['slope']:.4f} (bar {RATE_BAR}, "
         f"expected ≈ 3, predicted 2.6-3.6); pairwise rates "
         + ", ".join(f"{p:.4f}" for p in out["pairwise"]))
    _say(f"MAT21-3 [steady] max|Im theta| = {out['imag']:.3e} (predicted 0); "
         f"elapsed {time.perf_counter() - t0:.1f} s")
    return out


@pytest.fixture(scope="module")
def transient_p():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock3("transient")
    out = {"E": {}, "imag": 0.0}
    m = M_T
    msh, ft, ns, n, v_h, st = _rung_p(m, clk, DEGREE)
    out["imag"] = max(out["imag"], _max_imag(st.temperature))
    chi = _interior_chi(msh)
    _say(f"MAT21-3 [transient] rung m = {m} {tuple(ns)}: cells {n}; degree {DEGREE}; "
         f"t_C = {T_C:.3f} s, t* = {T_STAR:.3f} s")

    # (T) CN N = 40 with two BE half-steps, to t*.
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE)
    clk.lap(f"CN N={N_MAIN}")
    out["E"]["cn"] = res.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    th_t = res.temperature - T_INF
    ex_t, npts = _exact_p(msh, _trans_ex(T_STAR))
    clk.lap(f"transient series interpolation at t* ({npts} P3 points)")
    out["T"] = _rel_p(msh, th_t, ex_t, chi)
    _say(f"MAT21-3 (T) CN N={N_MAIN} at t*: interior rel L2 = {out['T']:.6e} "
         f"(bound {INTERIOR_BAND:.0e}, predicted 5e-5 to 9e-4; step 1 P1 m = 6 read 9.071528e-4); "
         f"whole {_rel_p(msh, th_t, ex_t):.6e} (printed); balance {res.max_balance_residual:.3e}")

    # (c2) Dirichlet z-roots in the transient comparand.
    ex_c2, _ = _exact_p(msh, _trans_ex(T_STAR, roots=(None, None, (np.arange(1, 41) - 0.5) * np.pi)))
    out["c2"] = _rel_p(msh, th_t, ex_c2, chi)
    _say(f"MAT21-3 control (c2) Dirichlet z-roots at t*: interior miss = {out['c2']:.4e} "
         f"(predicted 0.361, range 0.355-0.367; asserted > {INTERIOR_BAND:.0e})")
    clk.lap("control (c2)")

    # (S-slope) one BE step of DT_1 from theta = 0.
    res1 = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=DT_1, n_steps=1, theta=1.0, degree=DEGREE)
    clk.lap("BE one step DT_1")
    out["E"]["be_1"] = res1.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res1.temperature))
    th_bar = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=msh)))) / v_h
    ratio = RHO_C * th_bar / (DT_1 * Q)
    out["slope_dev"] = abs(ratio - 1.0)
    _say(f"MAT21-3 (S-slope) rho c theta_bar^1/(dt_1 q) = {ratio:.12f}; |ratio - 1| = "
         f"{out['slope_dev']:.4e} (bound {S2_BAND:.0e}, predicted 1.0144e-5 to 1.0146e-5); balance "
         f"{res1.max_balance_residual:.3e}")

    # (inf) backward Euler 20 x 5 t_C against the steady solve.
    res_inf = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=5.0 * T_C, n_steps=20, theta=1.0, degree=DEGREE)
    clk.lap("BE (inf) 20 x 5 t_C")
    out["E"]["be_inf"] = res_inf.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res_inf.temperature))
    # Step 1's argument form (test_mat21_box_fem.py:291).
    out["inf"] = _rel_p(msh, res_inf.temperature - T_INF, st.temperature - T_INF)
    _say(f"MAT21-3 (inf) BE 20 x 5 t_C vs steady: rel L2 = {out['inf']:.3e} "
         f"(bound {INF_BAND:.0e}, predicted ≈ 1e-13); balance {res_inf.max_balance_residual:.3e}")

    worst = max(out["E"], key=out["E"].get)
    _say("MAT21-3 (E) per run " + ", ".join(f"{k_}: {v:.3e}" for k_, v in out["E"].items())
         + f"; worst {out['E'][worst]:.3e} on {worst} (bound {E_BAND:.0e}, predicted <= 2e-11); "
         f"[transient] max|Im theta| = {out['imag']:.3e} (predicted 0); "
         f"elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL
    assert H == H_GEL
    assert RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs(steady_p):
    for r in steady_p["rungs"]:
        assert r["ndof"] == r["ndof_pred"], f"(Pdeg) m {r['m']}: {r['ndof']} != {r['ndof_pred']}"


@complex_only
def test_s1_p2_steady_rate(steady_p):
    assert steady_p["slope"] >= RATE_BAR, (
        f"(S1) slope {steady_p['slope']:.4f} < {RATE_BAR}; pairwise {steady_p['pairwise']}")


@complex_only
def test_s1_p2_steady_interior_finest(steady_p):
    assert steady_p["interior"] <= INTERIOR_BAND, f"(S1) interior {steady_p['interior']:.4e}"


@complex_only
def test_s2_p2_gauss_balance_every_rung(steady_p):
    for r in steady_p["rungs"]:
        assert r["gauss"] <= S2_GAUSS_BAND, f"(S2) m {r['m']}: {r['gauss']:.3e}"


@complex_only
def test_control_c1_p2_doubled_h_misses(steady_p):
    assert steady_p["c1"] > INTERIOR_BAND, f"(c1) {steady_p['c1']:.4e} inside band"


@complex_only
def test_control_cdeg_p1_worse_on_coarse_rung(steady_p):
    c = steady_p["cdeg"]
    assert c["e1_int"] > c["e2_int"], f"(c-deg) P1 {c['e1_int']:.4e} <= P2 {c['e2_int']:.4e}"


@complex_only
def test_t_p2_transient_interior(transient_p):
    assert transient_p["T"] <= INTERIOR_BAND, f"(T) interior {transient_p['T']:.4e}"


@complex_only
def test_slope_p2_calorimetric_first_step(transient_p):
    assert transient_p["slope_dev"] <= S2_BAND, f"(S-slope) {transient_p['slope_dev']:.4e}"


@complex_only
def test_e_p2_every_run_balances(transient_p):
    for name, r in transient_p["E"].items():
        assert r <= E_BAND, f"(E) {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
def test_inf_p2_long_backward_euler_is_steady(transient_p):
    assert transient_p["inf"] <= INF_BAND, f"(inf) {transient_p['inf']:.3e} > {INF_BAND:.0e}"


@complex_only
def test_control_c2_p2_dirichlet_z_roots_miss(transient_p):
    assert transient_p["c2"] > INTERIOR_BAND, f"(c2) {transient_p['c2']:.4e} inside band"
