"""Example (`EX-119`, `mat:20`): the F2182-sized gel box heated by an end-weighted source,
handed to the P2 thermal solvers three ways -- UFL, DG0 cell average, P2 nodal interpolant.

`MAT-21` step 8 (W2) gates P32's separable source q(x) = Q [1 - cos(pi x / A)] -- zero at the
box centre, 2 Q at the +-x ends, like end-ring SAR -- on step 3's m = 4 rung at degree 2.  Here
that rung is re-run with the gate's own calls (copied; its pytest fixtures cannot be called)
and three readings no example has are printed beside it: the fields in ParaView, the hot-end
profile along the x centre line, and the face-heat split under this source against the
uniform source on the same mesh.

Everything -- every band, every helper -- is imported from the gate modules
(``test_mat21_box_p2_source_fem`` and the modules it imports from).  Records live only in the
example-local ``REC`` dict, each commented with its ``log:line``; they are printed beside the
readings and never asserted by equality.

Asserted (imported bands only): (Pdeg) dofs; (S1-sep), (S1-dg0), (S1-p2s) steady interior;
(S2-p2s) Gauss; (T-p2s) CN interior at t*; (E-p2s) CN balance; (F0) six face tags partition the
Robin facets (integer counts); (F1) the six face powers sum to ``source_power`` within
``S2_GAUSS_BAND`` for the P2-source and the uniform runs.  Controls: (c-p2s) the P2-source run
against the uniform series misses (predicted 0.512), (cF) the +x face dropped misses
(predicted 0.04-0.25).
Not claimed: a coil or EM SAR source, whether a coil SAR can be P2-interpolated, the exact
series of the delivered DG0 source, any h-rate or other rung, F2182 compliance, perfusion, a
container wall, an external check.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:20 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_20_f2182_gel_box_end_weighted_source_combined.xdmf``
and ``..._x_centre_line.csv``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, la
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_mat21_box_p2_source_fem as g8  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    _cell_average_dg0_tag,
)
from tests.validation.test_mat21_box_closed_form import A_BOX, Q  # noqa: E402
from tests.validation.test_mat21_box_dg0_source_fem import STEP6_S1  # noqa: E402
from tests.validation.test_mat21_box_fem import (  # noqa: E402
    B_BOX,
    C_BOX,
    T_STAR,
    V_BOX,
    _interior_chi,
    _steady_ex,
)
from tests.validation.test_mat21_box_fem_p2 import (  # noqa: E402
    CMP_DEGREE,
    DEGREE,
    QDEG_P,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)
from tests.validation.test_mat21_box_separable_fem import (  # noqa: E402
    M4,
    PROFILES,
    SEP_TERMS,
    STEP3_M4,
    _point_vals,
    _sep_ex,
    _sep_floor,
    _sep_trans_ex,
)

INTERIOR_BAND = g8.INTERIOR_BAND
S2_GAUSS_BAND = g8.S2_GAUSS_BAND
E_BAND = g8.E_BAND
STEP7_S1 = g8.STEP7_S1
ROBIN_TAG = g8.ROBIN_TAG
RHO_C = g8.RHO_C
_allsum = g8._allsum
_max_imag = g8._max_imag

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_20_f2182_gel_box_end_weighted_source"

GATE = "20261016T033407Z_MAT-21-step8-W2-p2-source.log"

# Records, printed beside readings only; never asserted by equality.
REC = {
    "s1_p2s": 3.486949e-05,   # GATE :86
    "s2_p2s": 1.677e-14,      # GATE :85 (round-off floor; band assert only)
    "t_p2s": 1.388286e-04,    # GATE :92
    "e_p2s": 1.079e-13,       # GATE :93 (round-off floor; band assert only)
    "c_p2s": 5.1239e-01,      # GATE :87
    "a_p": 0.999999778,       # GATE :82
    "p2_minus_ufl": 1.168639e-07,   # GATE :94
    "dg0_minus_ufl": 3.390117e-04,  # GATE :94
    "th0_series": 2.741147,    # GATE :86
    "th09_series": 14.919189,  # GATE :86
}

A, B, C = A_BOX, B_BOX, C_BOX
LHALF = (A, B, C)
AXNAME = "xyz"
FACES = {1: (0, +1), 2: (0, -1), 3: (1, +1), 4: (1, -1), 5: (2, +1), 6: (2, -1)}


class _Clock:
    """Per-phase wall clock, printed on rank 0 with flush, prefix EX-119."""

    def __init__(self):
        self.t = time.perf_counter()
        self.t0 = self.t

    def lap(self, phase):
        now = time.perf_counter()
        if MPI.COMM_WORLD.rank == 0:
            print(f"EX-119 clock {phase}: {now - self.t:.1f} s", flush=True)
        self.t = now


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    def say(msg=""):
        if rank0:
            print(msg, flush=True)

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )
    assert comm.size == 2, f"comm.size {comm.size}"
    assert DEGREE == 2 and M4 == 4, (DEGREE, M4)

    say("=" * 72)
    say("EX-119 -- the F2182-sized gel box heated by an end-weighted source (UFL / DG0 / P2)")
    say("=" * 72)
    say(f"[records] beside {GATE}; DEGREE {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}, "
        f"m = {M4}, profiles {PROFILES}, SEP_TERMS {SEP_TERMS}, ranks {comm.size}")
    say(f"[bands] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, S2_GAUSS_BAND {S2_GAUSS_BAND:.0e}, "
        f"E_BAND {E_BAND:.0e}")
    clk = _Clock()
    imag = 0.0

    # ---- (i) uniform control (U), the gate's code path -------------------------------
    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    tdim = msh.topology.dim
    fdim = tdim - 1
    ds_r = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
    imap_c = msh.topology.index_map(tdim)
    n_cells_all = imap_c.size_local + imap_c.num_ghosts
    cells = np.arange(n_cells_all, dtype=np.int32)
    cell_tags = dmesh.meshtags(msh, tdim, cells, np.full(n_cells_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm)"},
        title="mat:20 -- F2182-sized gel box, Robin on all six faces, end-weighted source "
              "Q[1 - cos(pi x/A)]: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(0.0, -1.0, 0.0),
        comm=comm,
    )
    clk.lap("setup figure (no-op unless FEM_EM_SETUP_FIGURES=1)")

    ndof = {}
    ndof["u"] = int(st.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(st.temperature))
    th_u = st.temperature - T_INF
    ex_u, npts = _exact_p(msh, _steady_ex())
    chi = _interior_chi(msh)
    u_whole = _rel_p(msh, th_u, ex_u)
    u_int = _rel_p(msh, th_u, ex_u, chi)
    say(f"[anchor] (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {ndof['u']}; whole {u_whole:.6e} "
        f"(step 3 {STEP3_M4['whole']}); interior {u_int:.6e} (step 3 {STEP3_M4['interior']})")
    clk.lap("uniform series interpolation + readings")
    sp_u = float(np.real(st.source_power))
    conv_u = float(np.real(_allsum(comm, H_GEL * th_u * ds_r)))

    # ---- (ii) UFL-source separable steady solve ---------------------------------------
    xs = ufl.SpatialCoordinate(msh)
    q_sep = Q * (1.0 - ufl.cos(ufl.pi * xs[0] / A_BOX))
    ss = solve_steady_conduction(msh, ft, k=K_GEL, q=q_sep, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    ndof["s"] = int(ss.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(ss.temperature))
    th_s = ss.temperature - T_INF
    ex_s, npts_s = _exact_p(msh, _sep_ex())
    s6_int = _rel_p(msh, th_s, ex_s, chi)
    s6_whole = _rel_p(msh, th_s, ex_s)
    say(f"[anchor] (S1-sep) UFL-source steady interior {s6_int:.6e} <= {INTERIOR_BAND:.0e} "
        f"(step 6 {STEP6_S1['interior']}); whole {s6_whole:.6e} (step 6 {STEP6_S1['whole']})")
    clk.lap("UFL-source solve + readings")

    # ---- (iii) DG0 cell-average hand-off ----------------------------------------------
    n_loc = msh.topology.index_map(tdim).size_local
    ct = dmesh.meshtags(msh, tdim, np.arange(n_loc, dtype=np.int32),
                        np.ones(n_loc, dtype=np.int32))
    W, cbar = _cell_average_dg0_tag(msh, ct, 1, ufl.cos(ufl.pi * xs[0] / A_BOX))
    q_h = fem.Function(W, name="q_dg0")
    q_h.x.array[:] = Q * (1.0 - cbar)
    q_h.x.scatter_forward()
    sd = solve_steady_conduction(msh, ft, k=K_GEL, q=q_h, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    ndof["d"] = int(sd.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(sd.temperature))
    th_d = sd.temperature - T_INF
    d_int = _rel_p(msh, th_d, ex_s, chi)
    say(f"[anchor] (S1-dg0) DG0-source steady interior vs P32 series {d_int:.6e} <= "
        f"{INTERIOR_BAND:.0e} (step 7 {STEP7_S1['dg0_int']})")
    clk.lap("DG0-source solve + readings")

    # ---- (iv) P2 nodal interpolant of q ------------------------------------------------
    V2 = fem.functionspace(msh, ("Lagrange", DEGREE))
    q_p = fem.Function(V2, name="q_p2")
    q_p.interpolate(lambda x: Q * (1.0 - np.cos(np.pi * x[0] / A_BOX)))
    q_p.x.scatter_forward()
    dxq = ufl.dx(domain=msh, metadata=QDEG_P)
    a_p = float(np.real(_allsum(comm, (1.0 - q_p / Q) * ufl.cos(ufl.pi * xs[0] / A_BOX) * dxq))) \
        / (V_BOX / 2.0)
    say(f"[print] a_p = {a_p:.9f} (record {REC['a_p']}, {GATE}:82); 1 - a_p = {1.0 - a_p:+.4e}")

    sp2 = solve_steady_conduction(msh, ft, k=K_GEL, q=q_p, h=H_GEL, t_inf=T_INF,
                                  robin_tag=ROBIN_TAG, degree=DEGREE)
    ndof["p"] = int(sp2.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(sp2.temperature))
    th_p = sp2.temperature - T_INF
    conv_p = float(np.real(_allsum(comm, H_GEL * th_p * ds_r)))
    sp_p = float(np.real(sp2.source_power))
    gauss_p = abs(conv_p - sp_p) / sp_p
    p_int = _rel_p(msh, th_p, ex_s, chi)
    c_p2s = _rel_p(msh, th_p, ex_u, chi)
    pv = _point_vals(sp2.temperature, [[0.0, 0.0, 0.0], [0.9 * A_BOX, 0.0, 0.0]])
    say(f"[anchor] (S1-p2s) P2-source steady interior {p_int:.6e} <= {INTERIOR_BAND:.0e} "
        f"(record {REC['s1_p2s']:.6e}, {GATE}:86); theta_h(0) {pv[0]:.6f} K, theta_h(0.9A) "
        f"{pv[1]:.6f} K (series {REC['th0_series']:.6f} / {REC['th09_series']:.6f}, {GATE}:86)")
    say(f"[anchor] (S2-p2s) Gauss {gauss_p:.3e} <= {S2_GAUSS_BAND:.0e} (record {REC['s2_p2s']:.3e}, "
        f"{GATE}:85; round-off floor, band only); source_power/(Q V_h) - 1 = "
        f"{sp_p / (Q * v_h) - 1.0:+.3e}")
    say(f"[control] (c-p2s) P2-source run vs uniform series {c_p2s:.4e} > {INTERIOR_BAND:.0e} "
        f"(predicted 0.512, range [0.50, 0.53]; record {REC['c_p2s']:.4e}, {GATE}:87)")
    clk.lap("P2-source solve + readings")

    # ---- (v) CN N = N_MAIN to t* with the P2 source ------------------------------------
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_p, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE)
    imag = max(imag, _max_imag(res.temperature))
    e_bal = float(res.max_balance_residual)
    th_t = res.temperature - T_INF
    ex_t, _ = _exact_p(msh, _sep_trans_ex(T_STAR))
    t_int = _rel_p(msh, th_t, ex_t, chi)
    clk.lap(f"CN N={N_MAIN} (P2 source) + readings")
    say(f"[anchor] (T-p2s) CN N={N_MAIN} at t* = {T_STAR:.1f} s: interior {t_int:.6e} <= "
        f"{INTERIOR_BAND:.0e} (record {REC['t_p2s']:.6e}, {GATE}:92)")
    say(f"[anchor] (E-p2s) CN balance {e_bal:.3e} <= {E_BAND:.0e} (record {REC['e_p2s']:.3e}, "
        f"{GATE}:93; round-off floor, band only)")

    # ---- (vi) source-path differences, floor -------------------------------------------
    dp = _rel_p(msh, th_p, th_s, chi)
    dd = _rel_p(msh, th_d, th_s, chi)
    say(f"[print] source-path differences (interior, rel. to UFL run): P2 - UFL {dp:.6e} "
        f"(record {REC['p2_minus_ufl']:.6e}, {GATE}:94; predicted ≈ 1e-7); DG0 - UFL {dd:.6e} "
        f"(record {REC['dg0_minus_ufl']:.6e}, {GATE}:94; "
        f"predicted ≈ 3.4e-4)")
    say(f"[print] comparand floor 1 - prod S' at SEP_TERMS {_sep_floor():.4e}; max|Im theta| {imag:.3e}")

    pred = _dofs_pred(M4)
    say(f"[anchor] (Pdeg) P2 dofs uniform {ndof['u']} / UFL {ndof['s']} / DG0 {ndof['d']} / "
        f"P2-source {ndof['p']} (predicted {pred})")

    assert all(v == pred for v in ndof.values()), f"(Pdeg) {ndof} vs {pred}"
    assert s6_int <= INTERIOR_BAND, f"(S1-sep) {s6_int:.4e}"
    assert d_int <= INTERIOR_BAND, f"(S1-dg0) {d_int:.4e}"
    assert p_int <= INTERIOR_BAND, f"(S1-p2s) {p_int:.4e}"
    assert gauss_p <= S2_GAUSS_BAND, f"(S2-p2s) {gauss_p:.3e}"
    assert t_int <= INTERIOR_BAND, f"(T-p2s) {t_int:.4e}"
    assert e_bal <= E_BAND, f"(E-p2s) {e_bal:.3e}"
    assert c_p2s > INTERIOR_BAND, f"(c-p2s) did not miss: {c_p2s:.4e}"

    # ---- (F0) six face tags partition the Robin facets ---------------------------------
    msh.topology.create_connectivity(fdim, tdim)
    size_f = msh.topology.index_map(fdim).size_local
    idx_all, tag_all = [], []
    for tag, (d, sgn) in FACES.items():
        f = dmesh.locate_entities_boundary(
            msh, fdim,
            lambda x, d=d, sgn=sgn: np.isclose(x[d], sgn * LHALF[d], rtol=0, atol=1e-9))
        idx_all.append(np.asarray(f, dtype=np.int32))
        tag_all.append(np.full(f.size, tag, dtype=np.int32))
    idx = np.concatenate(idx_all)
    tg = np.concatenate(tag_all)
    order = np.argsort(idx, kind="stable")
    face_tags = dmesh.meshtags(msh, fdim, idx[order], tg[order])
    face_tags.name = "face_tags"
    nd = tuple(ns)  # (52, 32, 8); 2 * N_a * N_b triangles per face
    e1, e3, e5 = nd[1] * nd[2] * 2, nd[0] * nd[2] * 2, nd[0] * nd[1] * 2
    expected_counts = {1: e1, 2: e1, 3: e3, 4: e3, 5: e5, 6: e5}
    counts = {}
    for tag in FACES:
        loc = face_tags.find(tag)
        counts[tag] = comm.allreduce(int(np.count_nonzero(loc < size_f)), op=MPI.SUM)
    robin_loc = ft.find(ROBIN_TAG)
    n_robin = comm.allreduce(int(np.count_nonzero(robin_loc < size_f)), op=MPI.SUM)
    say()
    for tag, (d, sgn) in FACES.items():
        say(f"[F0] tag {tag} ({'+' if sgn > 0 else '-'}{AXNAME[d]}): owned facets {counts[tag]} "
            f"(expected {expected_counts[tag]})")
    say(f"[F0] sum of tagged facets {sum(counts.values())} vs Robin facets {n_robin} "
        f"(expected {sum(expected_counts.values())})")
    assert counts == expected_counts, f"(F0) per-tag counts {counts} vs {expected_counts}"
    assert sum(counts.values()) == sum(expected_counts.values()) == n_robin, (counts, n_robin)

    # ---- (F1) face powers, cosine and uniform sources ----------------------------------
    ds = ufl.Measure("ds", domain=msh, subdomain_data=face_tags)
    pw_p = {tag: float(np.real(_allsum(comm, H_GEL * th_p * ds(tag)))) for tag in FACES}
    pw_u = {tag: float(np.real(_allsum(comm, H_GEL * th_u * ds(tag)))) for tag in FACES}
    f1_p = abs(sum(pw_p.values()) / sp_p - 1.0)
    f1_u = abs(sum(pw_u.values()) / sp_u - 1.0)
    say(f"[F1] P2-source: |sum P_f / source_power - 1| = {f1_p:.3e} <= {S2_GAUSS_BAND:.0e}")
    say(f"[F1] uniform:   |sum P_f / source_power - 1| = {f1_u:.3e} <= {S2_GAUSS_BAND:.0e} "
        f"(Gauss reading {abs(conv_u - sp_u) / sp_u:.3e})")
    drop = sum(p for tag, p in pw_p.items() if tag != 1)
    cf = abs(drop / sp_p - 1.0)
    say(f"[control] (cF) +x face dropped: |sum_(f!=1) P_f / source_power - 1| = {cf:.5f} > "
        f"{S2_GAUSS_BAND:.0e} (predicted in [0.04, 0.25], not measured before; ceiling 0.5)")
    assert f1_p <= S2_GAUSS_BAND, f"(F1) P2-source {f1_p:.3e}"
    assert f1_u <= S2_GAUSS_BAND, f"(F1) uniform {f1_u:.3e}"
    assert cf > S2_GAUSS_BAND, f"(cF) did not miss: {cf:.3e}"
    say()
    say("[faces] pair | share of heat, cosine source | share, uniform source | change")
    for d, (t_pos, t_neg) in enumerate(((1, 2), (3, 4), (5, 6))):
        shp = (pw_p[t_pos] + pw_p[t_neg]) / sp_p
        shu = (pw_u[t_pos] + pw_u[t_neg]) / sp_u
        say(f"  ±{AXNAME[d]}  {shp:.6f} | {shu:.6f} | {shp - shu:+.6f}"
            + ("  (x-pair share predicted to rise; not asserted)" if d == 0 else ""))
    clk.lap("face tags and powers")

    # ---- x centre line, 129 points ------------------------------------------------------
    xl = np.linspace(-0.995 * A, 0.995 * A, 129)
    pts = np.zeros((xl.size, 3))
    pts[:, 0] = xl
    t_p = np.array(_point_vals(sp2.temperature, pts))
    t_d = np.array(_point_vals(sd.temperature, pts))
    t_uf = np.array(_point_vals(st.temperature, pts))
    t_se = np.real(_sep_ex()(pts[:, 0], pts[:, 1], pts[:, 2]))
    assert not np.any(np.isnan(t_p)) and not np.any(np.isnan(t_d)) and not np.any(np.isnan(t_uf))
    j = int(np.argmax(t_p))
    say()
    say(f"[line] x centre line, 129 points: max theta (P2 source) {t_p[j]:.6f} K at x = "
        f"{xl[j] / A:+.4f} A; theta(0) {pv[0]:.6f} K and theta(0.9A) {pv[1]:.6f} K beside series "
        f"{REC['th0_series']:.6f} / {REC['th09_series']:.6f} K ({GATE}:86)")
    k0 = int(np.argmin(np.abs(xl)))
    say(f"[line] theta at x = {xl[k0] / A:+.4f} A: P2 {t_p[k0]:.6f}, series {t_se[k0]:.6f}, DG0 "
        f"{t_d[k0]:.6f}, uniform-source FEM {t_uf[k0]:.6f} K")
    say(f"[line] max |P2 - series| along the line {np.max(np.abs(t_p - t_se)):.4e} K; "
        f"max |DG0 - P2| {np.max(np.abs(t_d - t_p)):.4e} K")
    clk.lap("centre line")

    # ---- outputs: CSV, then XDMF, then adopt_host_ownership -----------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_x_centre_line.csv"
        with open(csv_path, "w") as fh:
            fh.write("x_m,x_over_A,theta_p2src_K,theta_series_K,theta_dg0_K,theta_uniform_K\n")
            for i in range(xl.size):
                fh.write(f"{xl[i]:.9e},{xl[i] / A:.6f},{t_p[i]:.9e},{t_se[i]:.9e},"
                         f"{t_d[i]:.9e},{t_uf[i]:.9e}\n")
    comm.Barrier()

    V1 = fem.functionspace(msh, ("Lagrange", 1))

    def p1(name, src, offset=0.0):
        f = fem.Function(V1, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - offset
        f.x.scatter_forward()
        return f

    f_q = p1("q_p2", q_p)
    f_p = p1("theta_p2src_steady", sp2.temperature, T_INF)
    f_se = p1("theta_series_steady", ex_s)
    f_u = p1("theta_uniform_steady", st.temperature, T_INF)
    f_d = p1("theta_dg0_steady", sd.temperature, T_INF)
    f_t = p1("theta_tstar_p2src", res.temperature, T_INF)
    f_ts = p1("theta_tstar_series", ex_t)
    f_pd = fem.Function(V1, name="theta_p2src_minus_series")
    f_pd.x.array[:] = np.real(f_p.x.array) - np.real(f_se.x.array)
    f_dd = fem.Function(V1, name="theta_dg0_minus_p2src")
    f_dd.x.array[:] = np.real(f_d.x.array) - np.real(f_p.x.array)
    f_pd.x.scatter_forward()
    f_dd.x.scatter_forward()
    q_h.x.array[:] = np.real(q_h.x.array)
    fields = {}
    for f in (f_q, f_p, f_se, f_pd, f_dd, f_u, f_t, f_ts, q_h):
        fields[f.name] = f
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields, comm=comm,
        facet_tags=face_tags)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    clk.lap("outputs")
    if rank0:
        print(f"\n[output] {xdmf_path.name}: CG1 q_p2, theta_p2src_steady, theta_series_steady, "
              "theta_p2src_minus_series, theta_dg0_minus_p2src, theta_uniform_steady, "
              f"theta_tstar_p2src, theta_tstar_series (K, t* = {T_STAR:.1f} s), DG0 q_dg0, "
              "six face tags (theta_dg0_steady also written)")
        print(f"[output] {BASENAME}_x_centre_line.csv ({xl.size} rows)")
        print(f"[time] total {time.perf_counter() - clk.t0:.1f} s")
        print("[scope] one analytic source q = Q[1 - cos(pi x/A)] on the m = 4 rung at degree 2, -n 2; "
              "no coil or EM SAR, no h-rate, no F2182 compliance, no external check claimed")
        print("All assertions hold")


if __name__ == "__main__":
    main()
