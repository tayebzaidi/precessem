"""`MAT-21` step 20 FEM — the F2182 gel box with a meshed 6 mm acrylic wall
(§9 item 670, §10 P48b; brief `docs/planning/chunks/MAT-21.md` "§9 item 670").

No EM solve, no gmsh, no generator in `src`.  Mesh route (test-local): DolfinX
``create_box`` over the outer box with N_gel + 2 N_w cells per walled axis, then
every rank's geometry array (owned and ghost nodes) is remapped piecewise
linearly so grid planes fall exactly at ±A, ±B, ±C before any form is built.
Cells whose midpoint lies outside the gel box are wall (k_w, no source); k is a
DG0 field, q a DG0 indicator times Q.  Wall constants (d, k_w) imported from
item 660's module (nominal PMMA value, citation owed there).

W1 — the 1D-reduced geometry (the gate): the gel box in x, y (A, B ≫ C) with
the wall on ±z only, Robin on z = ±(C + d) only, the four side faces left out
of the Robin tag (natural Neumann), which is exactly item 660's 1D two-layer
slab.  Rungs m ∈ RUNGS: (13m, 8m, 2m + 2·N_w), N_w = m/2.  Asserted: (V) gel and
wall volumes = analytic (bar max(1e-12, 50·N·ε): the 2026-10-12 weekly's rule
(f) re-base of this mesh precondition, as `MAT-21` step 1); (S1) whole-domain
rel L2 LSQ slope ≥ ``RATE_BAR``; finest gel interior (|x|<0.9A, |y|<0.9B,
|z|<0.9C) ≤ ``INTERIOR_BAND``; (E) every rung's balance |∮h θ − qV_gel|/(qV_gel)
≤ ``E_BAND``.  Control (asserted miss): the finest FEM centre rise misses the
bare Robin slab series at h (no wall) by > ``INTERIOR_BAND``.

W2 — the six-wall box, one rung (print): centre rise of the meshed six-wall box
beside the bare gel box solved at h_eff on the same gel grid and the triple
series at h_eff and at h.  Asserted there only: the volume identity and the
balance (E); the centre figures are a print, no claim.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

import pytest

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction
from fem_em_solver.utils.analytical import (
    robin_box_steady_temperature,
    two_layer_slab_steady_temperature,
    wall_effective_heat_transfer,
)

from tests.validation.test_mat10_sar_to_heat_sphere import _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_finite_cylinder_fem import RATE_BAR, _rel
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    H,
    K,
    Q,
    SERIES_TERMS,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import T_INF
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, _bare_series

RUNGS = (2, 4, 6)
BASE = (13, 8, 2)
W2_RUNG = 4
ROBIN = 1
GEL, WALL = 1, 2
INTERP_CHUNK = 20000


def _n_wall(m):
    return m // 2


def _remap(s, L, d, n_gel, n_w):
    """Uniform grid on [−(L+d), L+d] (n_gel + 2 n_w cells) → planes at ±L."""
    if n_w == 0:
        return s
    tot = n_gel + 2 * n_w
    u = (s + (L + d)) / (2.0 * (L + d)) * tot  # grid-index coordinate
    out = np.where(
        u <= n_w, -(L + d) + u * (d / n_w),
        np.where(u <= n_w + n_gel, -L + (u - n_w) * (2.0 * L / n_gel),
                 L + (u - n_w - n_gel) * (d / n_w)))
    return out


def _walled_box(m, walls):
    """create_box + remap; ``walls`` is a per-axis bool.  Returns mesh, ns, nw."""
    comm = MPI.COMM_WORLD
    nw = _n_wall(m)
    ng = [b * m for b in BASE]
    ns = [ng[i] + (2 * nw if walls[i] else 0) for i in range(3)]
    outer = [HALF[i] + (D_WALL if walls[i] else 0.0) for i in range(3)]
    msh = dmesh.create_box(comm, [np.array([-o for o in outer]), np.array(outer)],
                           ns, cell_type=dmesh.CellType.tetrahedron)
    x = msh.geometry.x  # owned + ghost nodes on this rank
    for i in range(3):
        if walls[i]:
            x[:, i] = _remap(x[:, i].copy(), HALF[i], D_WALL, ng[i], nw)
    return msh, ns, nw


def _materials(msh):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    gel = ((np.abs(mids[:, 0]) < A_BOX) & (np.abs(mids[:, 1]) < B_BOX)
           & (np.abs(mids[:, 2]) < C_BOX))
    D0 = fem.functionspace(msh, ("DG", 0))
    dofs = D0.dofmap.list[:n, 0]
    chi = fem.Function(D0)
    chi.x.array[dofs] = gel.astype(float)
    kf = fem.Function(D0)
    kf.x.array[dofs] = np.where(gel, K, K_W)
    inter = fem.Function(D0)
    inter.x.array[dofs] = ((np.abs(mids[:, 0]) < 0.9 * A_BOX) & (np.abs(mids[:, 1]) < 0.9 * B_BOX)
                           & (np.abs(mids[:, 2]) < 0.9 * C_BOX)).astype(float)
    tags = np.where(gel, GEL, WALL).astype(np.int32)
    ct = dmesh.meshtags(msh, tdim, np.arange(n, dtype=np.int32), tags)
    return chi, kf, inter, ct


def _robin_facets(msh, outer, axes):
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)

    def mark(x):
        hit = np.zeros(x.shape[1], dtype=bool)
        for i in axes:
            hit |= np.isclose(np.abs(x[i]), outer[i], rtol=0.0, atol=1e-12)
        return hit

    f = np.sort(dmesh.locate_entities_boundary(msh, fdim, mark).astype(np.int32))
    return dmesh.meshtags(msh, fdim, f, np.full(f.size, ROBIN, dtype=np.int32))


def _volumes(msh, ct):
    comm = msh.comm
    dxc = ufl.Measure("dx", domain=msh, subdomain_data=ct)
    one = fem.Constant(msh, default_scalar_type(1.0))
    vg = float(np.real(_allsum(comm, one * dxc(GEL))))
    vw = float(np.real(_allsum(comm, one * dxc(WALL))))
    return vg, vw


def _centre(T):
    vals, ok = evaluate_vector_field_parallel(T, np.zeros((1, 3)))
    return float(np.real(np.ravel(vals[0])[0])) if bool(ok[0]) else float("nan")


def _slab_theta(z):
    return two_layer_slab_steady_temperature(
        np.asarray(z, dtype=float), half_thickness=C_BOX, wall=D_WALL, k=K, k_w=K_W,
        q=Q, h=H, t_inf=0.0)


def _exact_z(msh):
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    g = fem.Function(V2)
    g.interpolate(lambda x: np.real(_slab_theta(x[2])))
    return g


def _vol_bar(n):
    return max(1e-12, 50 * n * np.finfo(float).eps)


def _solve(msh, ft, kf, chi):
    return solve_steady_conduction(msh, ft, k=kf, q=Q * chi, h=H, t_inf=0.0,
                                   robin_tag=ROBIN, degree=1)


@pytest.fixture(scope="module")
def reduced():
    t0 = time.perf_counter()
    out = {"rungs": []}
    v_gel = 8.0 * A_BOX * B_BOX * C_BOX
    v_wall = 4.0 * A_BOX * B_BOX * 2.0 * D_WALL
    _say(f"\nMAT21-20F [W1] 1D-reduced: (A, B, C) = {HALF}, d = {D_WALL}, k = {K}, k_w = {K_W}, "
         f"h = {H}, q = {Q}; Robin on z = +-(C + d) only, side faces Neumann")
    for m in RUNGS:
        tl = time.perf_counter()
        msh, ns, nw = _walled_box(m, (False, False, True))
        n = msh.topology.index_map(3).size_global
        chi, kf, inter, ct = _materials(msh)
        ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
        vg, vw = _volumes(msh, ct)
        st = _solve(msh, ft, kf, chi)
        th = st.temperature
        ex = _exact_z(msh)
        e = _rel(msh, th, ex)
        ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN)
        conv = float(np.real(_allsum(msh.comm, H * th * ds)))
        bal = abs(conv - Q * v_gel) / (Q * v_gel)
        v_tot = vg + vw
        hm = (v_tot / n) ** (1.0 / 3.0)
        rec = {"m": m, "n": n, "h": hm, "e": e, "bal": bal,
               "dvg": abs(vg / v_gel - 1.0), "dvw": abs(vw / v_wall - 1.0)}
        _say(f"MAT21-20F [W1] m = {m} cells {tuple(ns)} N_w = {nw}: n = {n}; V_gel rel dev "
             f"{rec['dvg']:.3e}, V_wall rel dev {rec['dvw']:.3e} (bar {_vol_bar(n):.3e}); "
             f"h = {hm:.6e} m; whole rel L2 vs two-layer slab {e:.6e}; balance "
             f"|oint h theta - qV_gel|/(qV_gel) {bal:.3e} (solver: src {st.source_power:.9e} "
             f"conv {st.convected_power:.9e}); {time.perf_counter() - tl:.1f} s")
        out["rungs"].append(rec)
        if m == RUNGS[-1]:
            out["interior"] = _rel(msh, th, ex, inter)
            out["centre"] = _centre(th)
            out["centre_ex"] = float(_slab_theta(0.0))
    hs = np.array([r["h"] for r in out["rungs"]])
    es = np.array([r["e"] for r in out["rungs"]])
    out["slope"] = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
    pw = [float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1])) for i in range(len(hs) - 1)]
    _say(f"MAT21-20F [W1] (S1) LSQ slope {out['slope']:.4f} (bar {RATE_BAR}); pairwise "
         + ", ".join(f"{p:.4f}" for p in pw))
    _say(f"MAT21-20F [W1] (S1) finest interior rel L2 {out['interior']:.6e} (bound {INTERIOR_BAND:.0e})"
         f"; centre theta_h {out['centre']:.9f} K vs two-layer {out['centre_ex']:.9f} K "
         f"(rel {out['centre'] / out['centre_ex'] - 1.0:+.3e})")
    _say(f"MAT21-20F [W1] elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w1_volume_identity(reduced):
    for r in reduced["rungs"]:
        assert r["dvg"] <= _vol_bar(r["n"]) and r["dvw"] <= _vol_bar(r["n"]), r


def test_w1_steady_rate(reduced):
    assert reduced["slope"] >= RATE_BAR, reduced["slope"]


def test_w1_interior_finest(reduced):
    assert reduced["interior"] <= INTERIOR_BAND, reduced["interior"]


def test_w1_balance_every_rung(reduced):
    for r in reduced["rungs"]:
        assert r["bal"] <= E_BAND, r


def test_w1_control_bare_slab_at_h_misses(reduced):
    bare = float(_bare_series(0.0, H)[0]) - T_INF  # item 660's eigen-series at h, no wall
    miss = abs(reduced["centre"] - bare) / abs(bare)
    _say(f"MAT21-20F (control) FEM centre {reduced['centre']:.9f} K vs bare slab at h {bare:.9f} K: "
         f"rel miss {miss:.4e} (asserted > {INTERIOR_BAND:.0e})")
    assert miss > INTERIOR_BAND, miss


@pytest.fixture(scope="module")
def six_wall():
    t0 = time.perf_counter()
    m = W2_RUNG
    msh, ns, nw = _walled_box(m, (True, True, True))
    n = msh.topology.index_map(3).size_global
    chi, kf, _inter, ct = _materials(msh)
    outer = tuple(L + D_WALL for L in HALF)
    ft = _robin_facets(msh, outer, (0, 1, 2))
    vg, vw = _volumes(msh, ct)
    v_gel = 8.0 * A_BOX * B_BOX * C_BOX
    v_wall = 8.0 * np.prod(outer) - v_gel
    st = _solve(msh, ft, kf, chi)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN)
    conv = float(np.real(_allsum(msh.comm, H * st.temperature * ds)))
    bal = abs(conv - Q * v_gel) / (Q * v_gel)
    th_wall = _centre(st.temperature)
    # Bare gel box at h_eff on the same gel grid (no wall cells).
    h_eff = wall_effective_heat_transfer(H, D_WALL, K_W)
    ng = [b * m for b in BASE]
    mb = dmesh.create_box(MPI.COMM_WORLD, [-np.array(HALF), np.array(HALF)], ng,
                          cell_type=dmesh.CellType.tetrahedron)
    ftb = _robin_facets(mb, HALF, (0, 1, 2))
    stb = solve_steady_conduction(mb, ftb, k=K, q=Q, h=h_eff, t_inf=0.0, robin_tag=ROBIN, degree=1)
    th_bare = _centre(stb.temperature)
    kw = dict(half_widths=HALF, k=K, q=Q, t_inf=0.0, n_terms=SERIES_TERMS)
    ser_eff = float(np.real(robin_box_steady_temperature(0.0, 0.0, 0.0, h=h_eff, **kw)))
    ser_h = float(np.real(robin_box_steady_temperature(0.0, 0.0, 0.0, h=H, **kw)))
    out = {"n": n, "dvg": abs(vg / v_gel - 1.0), "dvw": abs(vw / v_wall - 1.0), "bal": bal}
    _say(f"\nMAT21-20F [W2] six-wall box m = {m} cells {tuple(ns)} N_w = {nw}: n = {n}; V_gel rel dev "
         f"{out['dvg']:.3e}, V_wall rel dev {out['dvw']:.3e} (bar {_vol_bar(n):.3e}); balance {bal:.3e}")
    _say(f"MAT21-20F [W2] (print) centre rise: six-wall meshed FEM {th_wall:.6f} K; bare gel box FEM at "
         f"h_eff {h_eff:.6f} (same gel grid {tuple(ng)}) {th_bare:.6f} K (h_eff-only miss "
         f"{th_bare / th_wall - 1.0:+.4e}); triple series at h_eff {ser_eff:.6f} K, at h {ser_h:.6f} K; "
         f"six-wall / series(h) {th_wall / ser_h:.6f}; elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w2_six_wall_volume_and_balance(six_wall):
    assert six_wall["dvg"] <= _vol_bar(six_wall["n"]), six_wall
    assert six_wall["dvw"] <= _vol_bar(six_wall["n"]), six_wall
    assert six_wall["bal"] <= E_BAND, six_wall
