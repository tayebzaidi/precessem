"""`MAT-15` step 2(b) — P1 on the phantom's own mesh against the finite Robin
cylinder series (§9 items 145 / 160, chain P2b; brief `docs/planning/chunks/MAT-15.md`).

No EM solve.  Mesh: ``_build(True, phantom_resolution=h_p)`` at `GEO-27`'s rungs
0.0075 / 0.005 / 0.00375 / 0.0025, tag-3 submesh, Robin on every exterior facet
(step 1's thermal domain).  Uniform q (2(a)'s ``Q``), gel k / h / ρc; comparand
2(a)'s ``robin_finite_cylinder_{steady,transient}_temperature`` interpolated into
P2 in numpy chunks, r = √(x² + y²).

Asserted: precondition (2(a)'s constants = step 1's; submesh centred/aligned);
(S1) steady whole-cylinder relative L2 LSQ slope against h_eff = (V_sub/N)^{1/3}
≥ 1.7, finest interior (r < 0.9a, |z| < 0.9L) ≤ ``INTERIOR_BAND``; (S2) Gauss
balance per rung ≤ 1e-9; (T) CN N = 40 with the Rannacher start to t* at 0.0025,
interior ≤ ``INTERIOR_BAND``; (E) every run's balance ≤ ``E_BAND``; (∞) BE
20 × 5 t_R = steady ≤ ``INF_BAND``.  Controls (rule (e), asserted only as
misses of ``INTERIOR_BAND``): (c1) comparand at h → 2h (steady), (c2) Dirichlet
radial roots in the transient comparand at t*.

Item 160 split the single module fixture into ``steady`` (four rungs: (S1),
(S2), precondition, (c1)) and ``transient`` (0.0025 rung alone: (T), (E), (∞),
(c2)), each phase clocked, so the two halves run as separate windows.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI
from scipy.special import jn_zeros

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_finite_cylinder_steady_temperature,
    robin_finite_cylinder_transient_temperature,
)
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _robin_tags, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, INF_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    T_R,
    T_STAR,
    A,
    _tag_submesh,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import KW
from tests.validation.test_mat15_finite_cylinder_closed_form import RHO_C as RHO_C_2A
from tests.validation.test_mat15_finite_cylinder_closed_form import A as A_2A
from tests.validation.test_mat15_finite_cylinder_closed_form import H as H_2A
from tests.validation.test_mat15_finite_cylinder_closed_form import K as K_2A
from tests.validation.test_mat15_finite_cylinder_closed_form import L, Q
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

RUNGS = (0.0075, 0.005, 0.00375, 0.0025)  # `GEO-27`'s phantom rungs
FINEST_CELLS = 58866
RATE_BAR = 1.7  # `MAT-9` step 1's bar
S2_GAUSS_BAND = 1e-9
# Centroid bar re-registered by measurement (2026-10-01 18:00 review ruling,
# `docs/planning/chunks/MAT-15.md` § item 145): the pre-registered 1e-9·a is
# unmeetable on an unstructured lateral faceting — measured |c| = 6.3e-8 m =
# 2.1e-6·a (`20261001T202818Z_MAT-15-step2b-W1.log:7163`).  1e-5·a is 4.8× that;
# a centroid offset δ moves the comparand by O(δ/a), four orders under
# ``INTERIOR_BAND``, and a 1 mm shift (3.3e-2·a) still fails it.
CENTROID_BAR = 1e-5
QDEG = {"quadrature_degree": 4}
STEADY_CHUNK = 20000  # (points × 400)
TRANSIENT_CHUNK = 1000  # (points × 200 × 200) in the einsum
V_CAD = np.pi * A**2 * 2.0 * L


def _steady_ex(**kw):
    def f(r, z):
        return robin_finite_cylinder_steady_temperature(r, z, t_inf=0.0, **{**KW, **kw})

    return f, STEADY_CHUNK


def _trans_ex(t, **kw):
    def f(r, z):
        return robin_finite_cylinder_transient_temperature(
            r, z, t, t_inf=0.0, rho_c=RHO_C, t_init=0.0, **{**KW, **kw}
        )

    return f, TRANSIENT_CHUNK


def _exact(sub, spec):
    fn, chunk = spec
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    g = fem.Function(V2)

    def ev(x):
        r = np.hypot(x[0], x[1])
        z = x[2]
        out = np.empty(r.shape)
        for s in range(0, r.size, chunk):
            out[s : s + chunk] = fn(r[s : s + chunk], z[s : s + chunk])
        return out

    g.interpolate(ev)
    return g


def _interior_chi(sub):
    tdim = sub.topology.dim
    imap = sub.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(sub, tdim, np.arange(n, dtype=np.int32))
    chi = fem.Function(fem.functionspace(sub, ("DG", 0)))
    inside = (np.hypot(mids[:, 0], mids[:, 1]) < 0.9 * A) & (np.abs(mids[:, 2]) < 0.9 * L)
    chi.x.array[:] = inside.astype(float)
    chi.x.scatter_forward()
    return chi


def _rel(sub, th_h, th_ex, chi=None):
    """sqrt(∫χ|θ_h − θ_ex|² / ∫χ|θ_ex|²), complex-safe."""
    dx = ufl.dx(domain=sub, metadata=QDEG)
    w = 1.0 if chi is None else chi
    d = th_h - th_ex
    e = _allsum(sub.comm, w * ufl.inner(d, d) * dx)
    n = _allsum(sub.comm, w * ufl.inner(th_ex, th_ex) * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


def _max_imag(T):
    n_own = T.function_space.dofmap.index_map.size_local
    loc = float(np.max(np.abs(np.imag(T.x.array[:n_own])))) if n_own else 0.0
    return T.function_space.mesh.comm.allreduce(loc, op=MPI.MAX)


def _geometry(sub, v_sub):
    comm = sub.comm
    x = sub.geometry.x
    zmin = comm.allreduce(float(x[:, 2].min()) if x.size else np.inf, op=MPI.MIN)
    zmax = comm.allreduce(float(x[:, 2].max()) if x.size else -np.inf, op=MPI.MAX)
    rmax = comm.allreduce(float(np.hypot(x[:, 0], x[:, 1]).max()) if x.size else 0.0, op=MPI.MAX)
    xs = ufl.SpatialCoordinate(sub)
    cen = np.array(
        [float(np.real(_allsum(comm, xs[i] * ufl.dx(domain=sub)))) / v_sub for i in range(3)]
    )
    return zmin, zmax, rmax, cen


class _Clock:
    """Per-phase wall clock, printed on rank 0 with flush (item 160)."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT15-2b [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _rung(hp, clk):
    """Build one rung's thermal submesh and solve the steady problem on it."""
    comm = MPI.COMM_WORLD
    msh, cell_tags, _ft, _d, t_build = _build(True, phantom_resolution=hp)
    sub, _pids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    n = sub.topology.index_map(3).size_global
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    clk.lap(f"build h_p={hp} (generator {t_build:.1f} s)")
    st = solve_steady_conduction(
        sub, ft, k=K_GEL, q=Q, h=H_GEL, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    clk.lap(f"steady solve h_p={hp}")
    return sub, ft, n, v_sub, st


@pytest.fixture(scope="module")
def steady():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock("steady")
    out = {"rungs": [], "imag": 0.0}
    _say(
        f"\nMAT15-2b Q = {Q} W/m^3, a = {A}, L = {L}, k = {K_GEL}, h = {H_GEL}, rho c = {RHO_C}; "
        f"t_R = {T_R:.3f} s, t* = {T_STAR:.3f} s"
    )
    for hp in RUNGS:
        sub, ft, n, v_sub, st = _rung(hp, clk)
        ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
        th_h = st.temperature - T_INF
        conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
        gauss = abs(conv - Q * v_sub) / (Q * v_sub)
        ex = _exact(sub, _steady_ex())
        clk.lap(f"steady series interpolation h_p={hp}")
        e_whole = _rel(sub, th_h, ex)
        h_eff = (v_sub / n) ** (1.0 / 3.0)
        out["imag"] = max(out["imag"], _max_imag(st.temperature))
        rec = {"hp": hp, "n": n, "v_sub": v_sub, "h_eff": h_eff, "e": e_whole, "gauss": gauss}
        _say(
            f"MAT15-2b rung h_p = {hp}: phantom cells {n}; V_sub = {v_sub:.12e} vs pi a^2 2L "
            f"{V_CAD:.12e} ({v_sub / V_CAD - 1.0:+.4e}, predicted ≈ -8.7e-4 at 0.0025); "
            f"h_eff = {h_eff:.5e} m; whole L2 rel = {e_whole:.6e}; (S2) Gauss "
            f"|oint h theta - qV|/qV = {gauss:.3e} (bound {S2_GAUSS_BAND:.0e})"
        )
        out["rungs"].append(rec)
        if hp != RUNGS[-1]:
            continue

        # Finest rung: precondition geometry, interior, centre, (c1).
        out["geom"] = _geometry(sub, v_sub)
        zmin, zmax, rmax, cen = out["geom"]
        _say(
            f"MAT15-2b geometry: z in [{zmin:.12e}, {zmax:.12e}] (L = {L}); max r = {rmax:.12e} "
            f"(a = {A}); centroid = ({cen[0]:.3e}, {cen[1]:.3e}, {cen[2]:.3e}), "
            f"|c|/a = {np.linalg.norm(cen) / A:.3e} (bar {CENTROID_BAR:.0e})"
        )
        chi = _interior_chi(sub)
        out["interior"] = _rel(sub, th_h, ex, chi)
        try:
            vals, ok = evaluate_vector_field_parallel(st.temperature, np.zeros((1, 3)))
            th_c = float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
        except Exception as exc:  # printed quantity only
            _say(f"MAT15-2b centre evaluation failed: {exc!r}")
            th_c = float("nan")
        th_c_ex = float(robin_finite_cylinder_steady_temperature(0.0, 0.0, t_inf=0.0, **KW))
        _say(
            f"MAT15-2b (S1) finest interior rel L2 = {out['interior']:.6e} (bound {INTERIOR_BAND:.0e}, "
            f"predicted 7.37e-4); theta_h(0) = {th_c:.6f} K vs series {th_c_ex:.6f} K "
            f"(P2a 2.983150; rel {th_c / th_c_ex - 1.0:+.3e})"
        )
        clk.lap("finest-rung interior + centre")
        ex_c1 = _exact(sub, _steady_ex(h=2.0 * H_GEL))
        out["c1"] = _rel(sub, th_h, ex_c1, chi)
        _say(
            f"MAT15-2b control (c1) comparand at 2h: steady interior miss = {out['c1']:.4e} "
            f"(predicted 0.739, asserted > {INTERIOR_BAND:.0e})"
        )
        clk.lap("control (c1)")

    hs = np.array([r["h_eff"] for r in out["rungs"]])
    es = np.array([r["e"] for r in out["rungs"]])
    out["slope"] = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
    out["pairwise"] = [
        float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1])) for i in range(len(hs) - 1)
    ]
    _say(
        f"MAT15-2b (S1) LSQ slope log e vs log h_eff = {out['slope']:.4f} (bar {RATE_BAR}, "
        f"predicted 2.15); pairwise rates " + ", ".join(f"{p:.4f}" for p in out["pairwise"])
    )
    _say(
        f"MAT15-2b [steady] max|Im theta| = {out['imag']:.3e} (predicted 0); "
        f"elapsed {time.perf_counter() - t0:.1f} s"
    )
    return out


def _run_transient(*, return_fields=False):
    """The `transient` fixture body, lifted to module level (rule (a), `EX-81`).

    Unset, byte-identical to the former fixture.  With ``return_fields=True`` the
    returned dict also carries ``sub``, ``ft``, ``v_sub``, ``st``, ``res``, ``ex_t``
    and ``chi`` (objects the body already built; additive only).
    """
    t0 = time.perf_counter()
    clk = _Clock("transient")
    out = {"E": {}, "imag": 0.0}
    hp = RUNGS[-1]
    sub, ft, n, v_sub, st = _rung(hp, clk)
    out["n"] = n
    out["imag"] = max(out["imag"], _max_imag(st.temperature))
    chi = _interior_chi(sub)
    _say(
        f"MAT15-2b [transient] rung h_p = {hp}: phantom cells {n}; "
        f"t_R = {T_R:.3f} s, t* = {T_STAR:.3f} s"
    )

    # (T) CN N = 40 with two BE half-steps, to t*.
    res = solve_transient_conduction(
        sub,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=Q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
    )
    clk.lap(f"CN N={N_MAIN}")
    out["E"]["cn"] = res.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    th_t = res.temperature - T_INF
    ex_t = _exact(sub, _trans_ex(T_STAR))
    clk.lap("transient series interpolation at t*")
    out["T"] = _rel(sub, th_t, ex_t, chi)
    _say(
        f"MAT15-2b (T) CN N={N_MAIN} at t*: interior rel L2 = {out['T']:.6e} "
        f"(bound {INTERIOR_BAND:.0e}, predicted 3.49e-4); whole {_rel(sub, th_t, ex_t):.6e} "
        f"(printed); balance {res.max_balance_residual:.3e}"
    )

    # (c2) Dirichlet radial roots in the transient comparand.
    ex_c2 = _exact(sub, _trans_ex(T_STAR, radial_roots=jn_zeros(0, 200)))
    out["c2"] = _rel(sub, th_t, ex_c2, chi)
    _say(
        f"MAT15-2b control (c2) Dirichlet radial roots at t*: interior miss = {out['c2']:.4e} "
        f"(predicted 0.779, asserted > {INTERIOR_BAND:.0e})"
    )
    clk.lap("control (c2)")

    # (∞) backward Euler 20 × 5 t_R against the steady solve.
    res_inf = solve_transient_conduction(
        sub,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=Q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=5.0 * T_R,
        n_steps=20,
        theta=1.0,
        degree=1,
    )
    clk.lap("BE (inf) 20 x 5 t_R")
    out["E"]["be_inf"] = res_inf.max_balance_residual
    out["imag"] = max(out["imag"], _max_imag(res_inf.temperature))
    # ‖T_be − T_st‖/‖θ_st‖ (MAT-12's :238 form). §9 item 168: the item-160 call passed
    # T_be − T_st as th_h, so it read 1 − O(agreement) (…034149Z…W3-transient.log:1882).
    out["inf"] = _rel(sub, res_inf.temperature - T_INF, st.temperature - T_INF)
    inf_ctl = _rel(sub, res_inf.temperature - st.temperature, st.temperature - T_INF)
    _say(
        f"MAT15-2b (inf) BE 20 x 5 t_R vs steady: rel L2 = {out['inf']:.3e} "
        f"(bound {INF_BAND:.0e}, predicted ≈ 1.9e-14); balance {res_inf.max_balance_residual:.3e}; "
        f"(inf-control) pre-fix expression = {inf_ctl:.16f} (printed, predicted 1.0)"
    )

    worst = max(out["E"], key=out["E"].get)
    _say(
        f"MAT15-2b (E) worst balance {out['E'][worst]:.3e} on {worst} (bound {E_BAND:.0e}); "
        f"[transient] max|Im theta| = {out['imag']:.3e} (predicted 0); "
        f"elapsed {time.perf_counter() - t0:.1f} s"
    )
    if return_fields:
        out.update(sub=sub, ft=ft, v_sub=v_sub, st=st, res=res, ex_t=ex_t, chi=chi)
    return out


@pytest.fixture(scope="module")
def transient():
    return _run_transient()


def test_precondition_2a_constants_are_step1s():
    assert A_2A == A
    assert L == PHANTOM_HEIGHT / 2.0
    assert K_2A == K_GEL
    assert H_2A == H_GEL
    assert RHO_C_2A == RHO_C


@complex_only
def test_precondition_submesh_centred_and_aligned(steady):
    zmin, zmax, rmax, cen = steady["geom"]
    assert abs(zmin + L) <= 1e-9 * L and abs(zmax - L) <= 1e-9 * L, (zmin, zmax)
    assert rmax <= A * (1.0 + 1e-9), rmax
    assert float(np.linalg.norm(cen)) <= CENTROID_BAR * A, cen


@complex_only
def test_finest_rung_cell_count(steady):
    assert steady["rungs"][-1]["n"] == FINEST_CELLS


@complex_only
def test_s1_steady_rate(steady):
    assert (
        steady["slope"] >= RATE_BAR
    ), f"(S1) slope {steady['slope']:.4f} < {RATE_BAR}; pairwise {steady['pairwise']}"


@complex_only
def test_s1_steady_interior_finest(steady):
    assert steady["interior"] <= INTERIOR_BAND, f"(S1) interior {steady['interior']:.4e} > 1 %"


@complex_only
def test_s2_gauss_balance_every_rung(steady):
    for r in steady["rungs"]:
        assert r["gauss"] <= S2_GAUSS_BAND, f"(S2) h_p {r['hp']}: {r['gauss']:.3e}"


@complex_only
def test_t_transient_interior(transient):
    assert transient["T"] <= INTERIOR_BAND, f"(T) interior {transient['T']:.4e} > 1 %"


@complex_only
def test_e_every_run_balances(transient):
    for name, r in transient["E"].items():
        assert r <= E_BAND, f"(E) {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
def test_inf_long_backward_euler_is_steady(transient):
    assert transient["inf"] <= INF_BAND, f"(inf) {transient['inf']:.3e} > {INF_BAND:.0e}"


@complex_only
def test_control_c1_doubled_h_misses(steady):
    assert steady["c1"] > INTERIOR_BAND, f"(c1) {steady['c1']:.4e} inside 1 %"


@complex_only
def test_control_c2_dirichlet_roots_miss(transient):
    assert transient["c2"] > INTERIOR_BAND, f"(c2) {transient['c2']:.4e} inside 1 %"
