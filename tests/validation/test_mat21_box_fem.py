"""`MAT-21` step 1 — P1 FEM on a structured ``create_box`` mesh of the F2182 gel
box against step 0's triple Robin eigen-series (§9 item 407; brief
`docs/planning/chunks/MAT-21.md`).

No EM solve.  Mesh: ``create_box`` over [−A, A] × [−B, B] × [−C, C] with
(13m, 8m, 2m) hexahedra split into 6 tetrahedra each, m ∈ {2, 3, 4, 6};
h_eff = (V/N)^{1/3} = 0.027/m m.  Robin on every exterior facet, uniform q,
gel k / h / ρc.  Comparand: step 0's ``robin_box_{steady,transient}_temperature``
(``n_terms=SERIES_TERMS``, ``t_inf=0``) interpolated into P2.

Asserted: (P0) step 0's constants = the FEM's; (P1) finest-rung bounds,
volume and centroid (volume bar max(1e-12, 50·N_cells·ε), the 2026-10-12
weekly's rule (f) re-base of a mesh precondition, not a physics band; §9 item
423 adds a planted 1 % volume-defect control that must miss it); (S1) steady whole-box relative-L2 LSQ slope ≥ ``RATE_BAR``,
finest interior ≤ ``INTERIOR_BAND``; (S2) Gauss balance per rung ≤
``S2_GAUSS_BAND``; (T) CN N = 40 at t* = 0.2 t_C interior ≤ ``INTERIOR_BAND``;
(S-slope) one BE step of Δt₁ = 1e-5 t_C, |ρc θ̄¹/(Δt₁ q) − 1| ≤ ``S2_BAND``;
(E) every run's balance ≤ ``E_BAND``; (∞) BE 20 × 5 t_C = steady ≤ ``INF_BAND``.
Controls (rule (e), asserted only as misses of ``INTERIOR_BAND``): (c1) the
comparand at h → 2h (steady, m = 6), (c2) Dirichlet z-roots at t*.
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
from fem_em_solver.utils.analytical import (
    robin_box_steady_temperature,
    robin_box_transient_temperature,
)

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
    _rel,
)
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    H,
    K,
    Q,
    SERIES_TERMS,
    T_C,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0

RUNGS = (2, 3, 4, 6)
BASE = (13, 8, 2)
T_STAR = 0.2 * T_C  # 2 821.5 s
DT_1 = 1.0e-5 * T_C  # 0.141 s, MAT-12's DT_1 = 1e-5 T_R pattern
V_BOX = 8.0 * A_BOX * B_BOX * C_BOX
INTERP_CHUNK = 2000  # points per series call (the function also chunks internally)


def _steady_ex(**kw):
    def f(x, y, z):
        return robin_box_steady_temperature(
            x, y, z, half_widths=HALF, k=K, q=Q, h=kw.get("h", H),
            t_inf=0.0, n_terms=SERIES_TERMS, roots=kw.get("roots"))
    return f


def _trans_ex(t, **kw):
    def f(x, y, z):
        return robin_box_transient_temperature(
            x, y, z, t, half_widths=HALF, k=K, q=Q, h=H, rho_c=RHO_C_0, t_inf=0.0,
            n_terms=SERIES_TERMS, roots=kw.get("roots"))
    return f


def _exact(msh, fn):
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    g = fem.Function(V2)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = np.real(
                fn(x[0, s:s + INTERP_CHUNK], x[1, s:s + INTERP_CHUNK], x[2, s:s + INTERP_CHUNK]))
        return out

    g.interpolate(ev)
    npts = msh.comm.allreduce(int(V2.dofmap.index_map.size_local), op=MPI.SUM)
    return g, npts


def _interior_chi(msh):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    chi = fem.Function(fem.functionspace(msh, ("DG", 0)))
    inside = ((np.abs(mids[:, 0]) < 0.9 * A_BOX) & (np.abs(mids[:, 1]) < 0.9 * B_BOX)
              & (np.abs(mids[:, 2]) < 0.9 * C_BOX))
    chi.x.array[:] = inside.astype(float)
    chi.x.scatter_forward()
    return chi


def _geometry(msh, v_h):
    comm = msh.comm
    x = msh.geometry.x
    lo = [comm.allreduce(float(x[:, i].min()) if x.size else np.inf, op=MPI.MIN) for i in range(3)]
    hi = [comm.allreduce(float(x[:, i].max()) if x.size else -np.inf, op=MPI.MAX) for i in range(3)]
    xs = ufl.SpatialCoordinate(msh)
    cen = np.array([float(np.real(_allsum(comm, xs[i] * ufl.dx(domain=msh)))) / v_h
                    for i in range(3)])
    return lo, hi, cen


class _Clock:
    """Per-phase wall clock, printed on rank 0 with flush (item 160)."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-1 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _rung(m, clk):
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
                                 robin_tag=ROBIN_TAG, degree=1)
    clk.lap(f"steady solve m={m}")
    return msh, ft, ns, n, v_h, st


@pytest.fixture(scope="module")
def steady():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock("steady")
    out = {"rungs": [], "imag": 0.0}
    _say(f"\nMAT21-1 Q = {Q} W/m^3, (A, B, C) = {HALF}, k = {K_GEL}, h = {H_GEL}, rho c = {RHO_C}; "
         f"t_C = {T_C:.3f} s, t* = {T_STAR:.3f} s, DT_1 = {DT_1:.5f} s")
    for m in RUNGS:
        msh, ft, ns, n, v_h, st = _rung(m, clk)
        ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
        th_h = st.temperature - T_INF
        conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
        gauss = abs(conv - Q * v_h) / (Q * v_h)
        ex, npts = _exact(msh, _steady_ex())
        clk.lap(f"steady series interpolation m={m} ({npts} P2 points)")
        e_whole = _rel(msh, th_h, ex)
        h_eff = (v_h / n) ** (1.0 / 3.0)
        out["imag"] = max(out["imag"], _max_imag(st.temperature))
        p1 = st.temperature.function_space.dofmap.index_map.size_global
        rec = {"m": m, "n": n, "v_h": v_h, "h_eff": h_eff, "e": e_whole, "gauss": gauss}
        _say(f"MAT21-1 rung m = {m} {tuple(ns)}: cells {n} (6*prod N = {6 * int(np.prod(ns))}); "
             f"P1 dofs {p1}; P2 points {npts}; V_h = {v_h:.15e} vs 8ABC {V_BOX:.15e} "
             f"({v_h / V_BOX - 1.0:+.3e}); h_eff = {h_eff:.6e} m (0.027/m = {0.027 / m:.6e}); "
             f"whole L2 rel = {e_whole:.6e}; (S2) Gauss |oint h theta - qV|/qV = {gauss:.3e} "
             f"(bound {S2_GAUSS_BAND:.0e}, predicted <= 1e-13)")
        out["rungs"].append(rec)
        if m != RUNGS[-1]:
            continue

        out["geom"] = _geometry(msh, v_h)
        out["v_h"] = v_h
        lo, hi, cen = out["geom"]
        _say(f"MAT21-1 (P1) geometry: lo = {lo}, hi = {hi}; centroid = ({cen[0]:.3e}, {cen[1]:.3e}, "
             f"{cen[2]:.3e}), |c|/A = {np.linalg.norm(cen) / A_BOX:.3e} (bar 1e-12)")
        chi = _interior_chi(msh)
        out["interior"] = _rel(msh, th_h, ex, chi)
        try:
            vals, ok = evaluate_vector_field_parallel(st.temperature, np.zeros((1, 3)))
            th_c = float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
        except Exception as exc:  # printed quantity only
            _say(f"MAT21-1 centre evaluation failed: {exc!r}")
            th_c = float("nan")
        th_c_ex = float(np.real(_steady_ex()(np.zeros(1), np.zeros(1), np.zeros(1))[0]))
        _say(f"MAT21-1 (S1) finest interior rel L2 = {out['interior']:.6e} (bound {INTERIOR_BAND:.0e}, "
             f"predicted 1e-3 to 8e-3); theta_h(0) = {th_c:.6f} K vs series {th_c_ex:.6f} K "
             f"(predicted 12.003947; rel {th_c / th_c_ex - 1.0:+.3e})")
        clk.lap("finest-rung interior + centre")
        ex_c1, _ = _exact(msh, _steady_ex(h=2.0 * H))
        out["c1"] = _rel(msh, th_h, ex_c1, chi)
        _say(f"MAT21-1 control (c1) comparand at 2h: steady interior miss = {out['c1']:.4e} "
             f"(predicted 0.60, range 0.55-0.65; asserted > {INTERIOR_BAND:.0e})")
        clk.lap("control (c1)")

    hs = np.array([r["h_eff"] for r in out["rungs"]])
    es = np.array([r["e"] for r in out["rungs"]])
    out["slope"] = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
    out["pairwise"] = [float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1]))
                       for i in range(len(hs) - 1)]
    _say(f"MAT21-1 (S1) whole L2 per rung " + ", ".join(f"{e:.4e}" for e in es)
         + " (predicted 2.5e-2 / 1.1e-2 / 6.2e-3 / 2.8e-3, x3)")
    _say(f"MAT21-1 (S1) LSQ slope log e vs log h_eff = {out['slope']:.4f} (bar {RATE_BAR}, "
         f"predicted 2.0); pairwise rates " + ", ".join(f"{p:.4f}" for p in out["pairwise"]))
    _say(f"MAT21-1 [steady] max|Im theta| = {out['imag']:.3e} (predicted 0); "
         f"elapsed {time.perf_counter() - t0:.1f} s")
    return out


@pytest.fixture(scope="module")
def transient():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock("transient")
    out = {"E": {}, "imag": 0.0}
    m = RUNGS[-1]
    msh, ft, ns, n, v_h, st = _rung(m, clk)
    out["imag"] = max(out["imag"], _max_imag(st.temperature))
    chi = _interior_chi(msh)
    _say(f"MAT21-1 [transient] rung m = {m} {tuple(ns)}: cells {n}; t_C = {T_C:.3f} s, "
         f"t* = {T_STAR:.3f} s")

    # (T) CN N = 40 with two BE half-steps, to t*.
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=1)
    clk.lap(f"CN N={N_MAIN}")
    out["E"]["cn"] = res.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    th_t = res.temperature - T_INF
    ex_t, npts = _exact(msh, _trans_ex(T_STAR))
    clk.lap(f"transient series interpolation at t* ({npts} P2 points)")
    out["T"] = _rel(msh, th_t, ex_t, chi)
    _say(f"MAT21-1 (T) CN N={N_MAIN} at t*: interior rel L2 = {out['T']:.6e} "
         f"(bound {INTERIOR_BAND:.0e}, predicted 5e-4 to 3e-3); whole {_rel(msh, th_t, ex_t):.6e} "
         f"(printed); balance {res.max_balance_residual:.3e}")

    # (c2) Dirichlet z-roots in the transient comparand.
    ex_c2, _ = _exact(msh, _trans_ex(T_STAR, roots=(None, None, (np.arange(1, 41) - 0.5) * np.pi)))
    out["c2"] = _rel(msh, th_t, ex_c2, chi)
    _say(f"MAT21-1 control (c2) Dirichlet z-roots at t*: interior miss = {out['c2']:.4e} "
         f"(predicted 0.365, range 0.33-0.40; asserted > {INTERIOR_BAND:.0e})")
    clk.lap("control (c2)")

    # (S-slope) one BE step of DT_1 from theta = 0.
    res1 = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=DT_1, n_steps=1, theta=1.0, degree=1)
    clk.lap("BE one step DT_1")
    out["E"]["be_1"] = res1.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res1.temperature))
    th_bar = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=msh)))) / v_h
    ratio = RHO_C * th_bar / (DT_1 * Q)
    out["slope_dev"] = abs(ratio - 1.0)
    _say(f"MAT21-1 (S-slope) rho c theta_bar^1/(dt_1 q) = {ratio:.12f}; |ratio - 1| = "
         f"{out['slope_dev']:.4e} (bound {S2_BAND:.0e}, predicted 1.0146e-5); balance "
         f"{res1.max_balance_residual:.3e}")

    # (inf) backward Euler 20 x 5 t_C against the steady solve.
    res_inf = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=5.0 * T_C, n_steps=20, theta=1.0, degree=1)
    clk.lap("BE (inf) 20 x 5 t_C")
    out["E"]["be_inf"] = res_inf.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res_inf.temperature))
    # The argument form of test_mat15_finite_cylinder_fem.py:297 (item 168).
    out["inf"] = _rel(msh, res_inf.temperature - T_INF, st.temperature - T_INF)
    inf_ctl = _rel(msh, res_inf.temperature - st.temperature, st.temperature - T_INF)
    _say(f"MAT21-1 (inf) BE 20 x 5 t_C vs steady: rel L2 = {out['inf']:.3e} "
         f"(bound {INF_BAND:.0e}, predicted ≈ 1e-13); balance {res_inf.max_balance_residual:.3e}; "
         f"(inf-control) pre-fix expression = {inf_ctl:.16f} (printed, predicted 1.0)")

    worst = max(out["E"], key=out["E"].get)
    _say(f"MAT21-1 (E) per run " + ", ".join(f"{k_}: {v:.3e}" for k_, v in out["E"].items())
         + f"; worst {out['E'][worst]:.3e} on {worst} (bound {E_BAND:.0e}, predicted <= 1e-12); "
         f"[transient] max|Im theta| = {out['imag']:.3e} (predicted 0); "
         f"elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL
    assert H == H_GEL
    assert RHO_C_0 == RHO_C


@complex_only
def test_precondition_box_geometry(steady):
    lo, hi, cen = steady["geom"]
    for i, L in enumerate(HALF):
        assert abs(lo[i] + L) <= 1e-12 * L and abs(hi[i] - L) <= 1e-12 * L, (i, lo[i], hi[i])
    # Volume bar re-based by the 2026-10-12 weekly's rule (f) ruling
    # (PROJECT_PLAN.md:4575; §9 item 423): item 407 measured 1.659e-12 at
    # N = 269 568 (W1 log 20261013T004905Z, :111) against the old 1e-12 —
    # summation round-off over N cells, not a geometry defect.  The bar is a
    # mesh precondition, not a physics band; bounds and centroid keep 1e-12.
    n_cells = steady["rungs"][-1]["n"]  # index_map(3).size_global
    vol_bar = max(1e-12, 50 * n_cells * np.finfo(float).eps)
    vol_dev = abs(steady["v_h"] / V_BOX - 1.0)
    _say(f"MAT21-1 (P1) volume |V_h/8ABC - 1| = {vol_dev:.3e} (bar max(1e-12, 50 N eps) = "
         f"{vol_bar:.3e} at N = {n_cells}; rule (f) re-base, PROJECT_PLAN.md:4575)")
    assert vol_dev <= vol_bar, (steady["v_h"], vol_bar)
    assert float(np.linalg.norm(cen)) <= 1e-12 * A_BOX, cen


@complex_only
def test_control_p1_planted_volume_defect_misses():
    """Planted 1 % volume defect (x half-width × 1.01) must miss the re-based bar."""
    comm = MPI.COMM_WORLD
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(comm, [np.array([-1.01 * A_BOX, -B_BOX, -C_BOX]),
                                  np.array([1.01 * A_BOX, B_BOX, C_BOX])],
                           ns, cell_type=dmesh.CellType.tetrahedron)
    one = fem.Constant(msh, default_scalar_type(1.0))
    v_planted = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    n = msh.topology.index_map(3).size_global
    bar = max(1e-12, 50 * n * np.finfo(float).eps)
    miss = abs(v_planted / V_BOX - 1.0)
    _say(f"MAT21-1 control (P1-planted) 1.01 A box: |V/8ABC - 1| = {miss:.6e} "
         f"(predicted 1.000e-2; asserted > bar); bar = {bar:.4e} at n = {n} "
         f"(predicted 2.993e-9 at n = 269568); separation {miss / bar:.3e}x "
         f"(predicted ≈ 3.3e6x, printed only)")
    assert miss > bar, (miss, bar)


@complex_only
def test_s1_steady_rate(steady):
    assert steady["slope"] >= RATE_BAR, (
        f"(S1) slope {steady['slope']:.4f} < {RATE_BAR}; pairwise {steady['pairwise']}")


@complex_only
def test_s1_steady_interior_finest(steady):
    assert steady["interior"] <= INTERIOR_BAND, f"(S1) interior {steady['interior']:.4e} > 1 %"


@complex_only
def test_s2_gauss_balance_every_rung(steady):
    for r in steady["rungs"]:
        assert r["gauss"] <= S2_GAUSS_BAND, f"(S2) m {r['m']}: {r['gauss']:.3e}"


@complex_only
def test_control_c1_doubled_h_misses(steady):
    assert steady["c1"] > INTERIOR_BAND, f"(c1) {steady['c1']:.4e} inside 1 %"


@complex_only
def test_t_transient_interior(transient):
    assert transient["T"] <= INTERIOR_BAND, f"(T) interior {transient['T']:.4e} > 1 %"


@complex_only
def test_slope_calorimetric_first_step(transient):
    assert transient["slope_dev"] <= S2_BAND, f"(S-slope) {transient['slope_dev']:.4e} > {S2_BAND:.0e}"


@complex_only
def test_e_every_run_balances(transient):
    for name, r in transient["E"].items():
        assert r <= E_BAND, f"(E) {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
def test_inf_long_backward_euler_is_steady(transient):
    assert transient["inf"] <= INF_BAND, f"(inf) {transient['inf']:.3e} > {INF_BAND:.0e}"


@complex_only
def test_control_c2_dirichlet_z_roots_miss(transient):
    assert transient["c2"] > INTERIOR_BAND, f"(c2) {transient['c2']:.4e} inside 1 %"
