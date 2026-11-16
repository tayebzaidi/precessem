"""Example (`EX-129`, `mat:23`): the F2182-sized gel box heated by the transverse eddy-SAR shape
q = 3Q (y^2 + z^2) / (B^2 + C^2) -- hot side walls, heat leaving through +-y.

`MAT-21` step 9 (W2) gates this long-prism circularly-polarised-B1 SAR shape on step 3's m = 4
rung at degree 2, steady and CN to t*.  The weights (w_y = 2.868293, w_z = 0.131707) are physics,
not tuning.  Here the rung is re-run with the gate's own calls (copied; its pytest fixture cannot
be called) and four readings no example has are printed beside it: a source that varies along
the two transverse axes, the y centre line, the face-heat split moving toward +-y, and the thin-
axis (z^2) term as a field in ParaView.

Everything -- every band, every helper -- is imported from the gate modules
(``tests.validation.test_mat21_box_quadratic_fem`` as ``g9`` and the modules it imports from).
Records live only in the example-local ``REC`` dict, each commented with its ``log:line``; they
are printed beside the readings and never asserted by equality.

Asserted (imported bands only): (Pdeg) dofs; (S1-quad) steady interior; (S2-quad) Gauss;
(T-quad) CN interior at t*; (E-quad) CN balance; (U) the uniform control reproduces step 3's
digits; (F0) six face tags partition the Robin facets (integer counts); (F1) the six face
powers sum to ``source_power`` within ``S2_GAUSS_BAND`` for both runs.  Backed controls
(each asserted > INTERIOR_BAND): (c-src), (c-src'), (c-z), (c-src-T).  Example-side control (cF):
the +y face dropped misses.
Not claimed: the box's own EM SAR (end-face current return near x = +-A, wavelength effects at
64 / 128 MHz), no EM is solved and no coil is involved, linear polarisation or any other weight
pair, a P2 / DG0 hand-off, the product with an axial cosine, any h-rate, another rung or degree 1,
an absolute SAR or temperature for a real scanner, F2182 compliance or probe placement advice,
perfusion, a container wall, an external check.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:23 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_23_f2182_gel_box_transverse_eddy_source_combined.xdmf``
and ``..._y_centre_line.csv``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_mat21_box_quadratic_fem as g9  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
)
from tests.validation.test_mat21_box_closed_form import (  # noqa: E402
    A_BOX,
    Q,
    SERIES_TERMS,
)
from tests.validation.test_mat21_box_fem import (  # noqa: E402
    B_BOX,
    C_BOX,
    T_STAR,
    V_BOX,
    _interior_chi,
    _steady_ex,
    _trans_ex,
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
from tests.validation.test_mat21_box_quadratic_closed_form import (  # noqa: E402
    QUAD_TERMS,
    WEIGHTS,
    WEIGHTS_NOZ,
)
from tests.validation.test_mat21_box_separable_fem import M4, STEP3_M4, _point_vals  # noqa: E402

INTERIOR_BAND = g9.INTERIOR_BAND
S2_GAUSS_BAND = g9.S2_GAUSS_BAND
E_BAND = g9.E_BAND
ROBIN_TAG = g9.ROBIN_TAG
RHO_C = g9.RHO_C
_allsum = g9._allsum
_max_imag = g9._max_imag

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_23_f2182_gel_box_transverse_eddy_source"

GATE = "20261016T171904Z_MAT-21-step9-W2-quadratic-fem.log"
W1 = "20261016T171632Z_MAT-21-step9-W1-closed-form.log"
EX119 = "20261016T112833Z_EX-119.log"

# Records, printed beside readings only; never asserted by equality.
REC = {
    "u_whole": 2.489307e-05, "u_int": 1.669359e-05,            # GATE :70
    "s1_int": 7.219878e-05, "s1_whole": 9.086572e-05,           # GATE :76
    "th0": 4.276604, "th09y": 13.337975, "th09x": 2.682175,     # GATE :76
    "th0_se": 4.276700, "th09y_se": 13.337014, "th09x_se": 2.682231,  # GATE :76
    "s2": 9.775e-14,                                            # GATE :73 (round-off floor)
    "sp_rel": -1.350e-13,                                       # GATE :73
    "c_src": 4.4856e-01, "c_srcp": 5.1471e-01, "c_z": 4.7893e-02,  # GATE :77-79
    "t_int": 2.694340e-04, "t_whole": 4.310154e-04,             # GATE :83
    "tv0": 0.063442, "tv09y": 2.893351,                         # GATE :83
    "tv0_se": 0.063466, "tv09y_se": 2.892661,                   # GATE :83
    "e": 2.285e-13,                                             # GATE :84 (round-off floor)
    "c_srcT": 6.9562e-01,                                       # GATE :85
    "floor": 1.2484e-06,                                        # GATE :86
    "q1": 2.904e-15, "q1b": 2.492e-15,                          # W1 :236
    "q5y": 2.484e-09, "q5z": 1.460e-11,                         # W1 :254
    "uni_y_share": 0.126246,                                    # EX119 :81
}

A, B, C = A_BOX, B_BOX, C_BOX
LHALF = (A, B, C)
AXNAME = "xyz"
FACES = {1: (0, +1), 2: (0, -1), 3: (1, +1), 4: (1, -1), 5: (2, +1), 6: (2, -1)}


class _Clock:
    """Per-phase wall clock, printed on rank 0 with flush, prefix EX-129."""

    def __init__(self):
        self.t = time.perf_counter()
        self.t0 = self.t

    def lap(self, phase):
        now = time.perf_counter()
        if MPI.COMM_WORLD.rank == 0:
            print(f"EX-129 clock {phase}: {now - self.t:.1f} s", flush=True)
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
    say("EX-129 -- the F2182-sized gel box under the transverse eddy source 3Q(y^2 + z^2)/(B^2 + C^2)")
    say("=" * 72)
    say(f"[records] beside {GATE}; DEGREE {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}, "
        f"m = {M4}, weights {WEIGHTS}, QUAD_TERMS {QUAD_TERMS}, SERIES_TERMS {SERIES_TERMS}, "
        f"ranks {comm.size}")
    say(f"[bands] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, S2_GAUSS_BAND {S2_GAUSS_BAND:.0e}, "
        f"E_BAND {E_BAND:.0e}")
    say(f"[provenance] closed form checked by quadrature and two reduced routes in {W1}: (Q1) "
        f"{REC['q1']:.3e} / (Q1b) {REC['q1b']:.3e} (:236); (Q5) y-route {REC['q5y']:.3e} / z-route "
        f"{REC['q5z']:.3e} (:254); the numpy module is not re-run here")
    clk = _Clock()
    imag = 0.0

    # ---- (i) uniform control (U), the gate's code path (g9:125-134) -----------------------
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
        title="mat:23 -- F2182-sized gel box, Robin on all six faces, transverse source "
              "3Q(y^2 + z^2)/(B^2 + C^2): x = 0 slice",
        slice_normal=(1.0, 0.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(-1.0, 0.0, 0.0),
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
        f"(step 3 {STEP3_M4['whole']}, {GATE}:70 record {REC['u_whole']:.6e}); interior {u_int:.6e} "
        f"(step 3 {STEP3_M4['interior']}, record {REC['u_int']:.6e})")
    clk.lap("uniform series interpolation + readings")
    sp_u = float(np.real(st.source_power))
    conv_u = float(np.real(_allsum(comm, H_GEL * th_u * ds_r)))

    # ---- (ii) quadratic steady solve (g9:141-144) -----------------------------------------
    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    sq = solve_steady_conduction(msh, ft, k=K_GEL, q=q_quad, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("quadratic steady solve")
    ndof["q"] = int(sq.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(sq.temperature))
    th_q = sq.temperature - T_INF
    conv = float(np.real(_allsum(comm, H_GEL * th_q * ds_r)))
    sp = float(np.real(sq.source_power))
    gauss = abs(conv - sp) / sp
    say(f"[anchor] (S2-quad) Gauss {gauss:.3e} <= {S2_GAUSS_BAND:.0e} (record {REC['s2']:.3e}, "
        f"{GATE}:73; round-off floor, band only); source_power/(Q V_h) - 1 = "
        f"{sp / (Q * v_h) - 1.0:+.3e} (record {REC['sp_rel']:+.3e}, {GATE}:73)")
    ex_q, npts_q = _exact_p(msh, g9._quad_ex())
    clk.lap(f"quadratic series interpolation ({npts_q} P3 points, {QUAD_TERMS})")
    q_int = _rel_p(msh, th_q, ex_q, chi)
    q_whole = _rel_p(msh, th_q, ex_q)
    ex_z, _ = _exact_p(msh, g9._quad_ex(WEIGHTS_NOZ, SERIES_TERMS))
    clk.lap(f"no-z series interpolation ({SERIES_TERMS})")
    c_src = _rel_p(msh, th_q, ex_u, chi)
    c_srcp = _rel_p(msh, th_u, ex_q, chi)
    c_z = _rel_p(msh, th_q, ex_z, chi)
    pv = _point_vals(sq.temperature, [[0.0, 0.0, 0.0], [0.0, 0.9 * B_BOX, 0.0],
                                      [0.9 * A_BOX, 0.0, 0.0]])
    say(f"[anchor] (S1-quad) steady interior vs quadratic series {q_int:.6e} <= {INTERIOR_BAND:.0e} "
        f"(record {REC['s1_int']:.6e}, {GATE}:76); whole {q_whole:.6e} (record "
        f"{REC['s1_whole']:.6e})")
    say(f"[print] steady points: theta_h(0) {pv[0]:.6f} K (record {REC['th0']:.6f}, series "
        f"{REC['th0_se']:.6f}), theta_h(0,0.9B,0) {pv[1]:.6f} K (record {REC['th09y']:.6f}, series "
        f"{REC['th09y_se']:.6f}), theta_h(0.9A,0,0) {pv[2]:.6f} K (record {REC['th09x']:.6f}, "
        f"series {REC['th09x_se']:.6f}) [{GATE}:76]")
    say(f"[control] (c-src) quadratic FEM vs uniform series {c_src:.4e} > {INTERIOR_BAND:.0e} "
        f"(backed, record {REC['c_src']:.4e}, {GATE}:77)")
    say(f"[control] (c-src') uniform FEM vs quadratic series {c_srcp:.4e} > {INTERIOR_BAND:.0e} "
        f"(backed, record {REC['c_srcp']:.4e}, {GATE}:78)")
    say(f"[control] (c-z) quadratic FEM vs no-z series {c_z:.4e} > {INTERIOR_BAND:.0e} "
        f"(backed, record {REC['c_z']:.4e}, {GATE}:79; ratio to band {c_z / INTERIOR_BAND:.2f}; "
        f"no factor bar posed)")
    clk.lap("quadratic steady readings")

    # ---- (iii) CN N = N_MAIN to t* (g9:180-183) --------------------------------------------
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_quad, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE)
    clk.lap(f"CN N={N_MAIN} (quadratic source)")
    ndof["t"] = int(res.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(res.temperature))
    e_bal = float(res.max_balance_residual)
    th_t = res.temperature - T_INF
    ex_t, _ = _exact_p(msh, g9._quad_trans_ex(T_STAR))
    ex_ut, _ = _exact_p(msh, _trans_ex(T_STAR))
    t_int = _rel_p(msh, th_t, ex_t, chi)
    t_whole = _rel_p(msh, th_t, ex_t)
    c_srcT = _rel_p(msh, th_t, ex_ut, chi)
    clk.lap("transient interpolation + readings")
    tv = _point_vals(res.temperature, [[0.0, 0.0, 0.0], [0.0, 0.9 * B_BOX, 0.0]])
    say(f"[anchor] (T-quad) CN N={N_MAIN} at t* = {T_STAR:.1f} s: interior {t_int:.6e} <= "
        f"{INTERIOR_BAND:.0e} (record {REC['t_int']:.6e}, {GATE}:83); whole {t_whole:.6e} "
        f"(record {REC['t_whole']:.6e})")
    say(f"[print] transient points: theta_h(0, t*) {tv[0]:.6f} K (record {REC['tv0']:.6f}, series "
        f"{REC['tv0_se']:.6f}), theta_h(0,0.9B,0, t*) {tv[1]:.6f} K (record {REC['tv09y']:.6f}, "
        f"series {REC['tv09y_se']:.6f}) [{GATE}:83]")
    say(f"[anchor] (E-quad) CN balance {e_bal:.3e} <= {E_BAND:.0e} (record {REC['e']:.3e}, "
        f"{GATE}:84; round-off floor, band only)")
    say(f"[control] (c-src-T) CN quadratic vs uniform transient series {c_srcT:.4e} > "
        f"{INTERIOR_BAND:.0e} (backed, record {REC['c_srcT']:.4e}, {GATE}:85)")

    fl = g9._quad_floor()
    say(f"[print] comparand floor 1 - target at QUAD_TERMS {fl:.4e} (record {REC['floor']:.4e}, "
        f"{GATE}:86); (S1-quad)/floor {q_int / fl:.1f} (record 57.8); max|Im theta| {imag:.3e}; "
        f"V_BOX/V_h - 1 {V_BOX / v_h - 1.0:+.3e}")

    pred = _dofs_pred(M4)
    say(f"[anchor] (Pdeg) P2 dofs uniform {ndof['u']} / quadratic {ndof['q']} / CN {ndof['t']} "
        f"(predicted {pred})")

    assert pred == 116025, pred
    assert all(v == pred for v in ndof.values()), f"(Pdeg) {ndof} vs {pred}"
    assert f"{u_whole:.6e}" == STEP3_M4["whole"], f"(U) whole {u_whole:.6e}"
    assert f"{u_int:.6e}" == STEP3_M4["interior"], f"(U) interior {u_int:.6e}"
    assert q_int <= INTERIOR_BAND, f"(S1-quad) {q_int:.4e}"
    assert gauss <= S2_GAUSS_BAND, f"(S2-quad) {gauss:.3e}"
    assert t_int <= INTERIOR_BAND, f"(T-quad) {t_int:.4e}"
    assert e_bal <= E_BAND, f"(E-quad) {e_bal:.3e}"
    assert c_src > INTERIOR_BAND, f"(c-src) did not miss: {c_src:.4e}"
    assert c_srcp > INTERIOR_BAND, f"(c-src') did not miss: {c_srcp:.4e}"
    assert c_z > INTERIOR_BAND, f"(c-z) did not miss: {c_z:.4e}"
    assert c_srcT > INTERIOR_BAND, f"(c-src-T) did not miss: {c_srcT:.4e}"

    # ---- (F0) six face tags partition the Robin facets (20_...py:307-338) ------------------
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
        f"(expected {sum(expected_counts.values())}; record 512/512/832/832/3328/3328 = 9344, "
        f"{EX119}:68-74)")
    assert counts == expected_counts, f"(F0) per-tag counts {counts} vs {expected_counts}"
    assert counts == {1: 512, 2: 512, 3: 832, 4: 832, 5: 3328, 6: 3328}, f"(F0) literal {counts}"
    assert sum(counts.values()) == sum(expected_counts.values()) == n_robin == 9344, (counts, n_robin)

    # ---- (F1) face powers, quadratic and uniform sources (20_...py:340-355) -----------------
    ds = ufl.Measure("ds", domain=msh, subdomain_data=face_tags)
    pw_q = {tag: float(np.real(_allsum(comm, H_GEL * th_q * ds(tag)))) for tag in FACES}
    pw_u = {tag: float(np.real(_allsum(comm, H_GEL * th_u * ds(tag)))) for tag in FACES}
    f1_q = abs(sum(pw_q.values()) / sp - 1.0)
    f1_u = abs(sum(pw_u.values()) / sp_u - 1.0)
    say(f"[F1] quadratic: |sum P_f / source_power - 1| = {f1_q:.3e} <= {S2_GAUSS_BAND:.0e}")
    say(f"[F1] uniform:   |sum P_f / source_power - 1| = {f1_u:.3e} <= {S2_GAUSS_BAND:.0e} "
        f"(Gauss reading {abs(conv_u - sp_u) / sp_u:.3e})")
    drop_q = sum(p for tag, p in pw_q.items() if tag != 3)
    drop_u = sum(p for tag, p in pw_u.items() if tag != 3)
    cf = abs(drop_q / sp - 1.0)
    cf_u = abs(drop_u / sp_u - 1.0)
    say(f"[control] (cF) +y face dropped: |sum_(f!=3) P_f / source_power - 1| = {cf:.5f} > "
        f"{S2_GAUSS_BAND:.0e} (predicted in [0.075, 0.20], not measured before; ceiling 0.5; "
        f"in range: {0.075 <= cf <= 0.20}); same miss on the uniform run {cf_u:.6f} "
        f"(half of {EX119}:81 = {REC['uni_y_share'] / 2.0:.6f})")
    assert f1_q <= S2_GAUSS_BAND, f"(F1) quadratic {f1_q:.3e}"
    assert f1_u <= S2_GAUSS_BAND, f"(F1) uniform {f1_u:.3e}"
    assert cf > S2_GAUSS_BAND, f"(cF) did not miss: {cf:.3e}"
    assert cf < 0.5, f"(cF) above ceiling: {cf:.3e}"
    say()
    say("[faces] pair | share of heat, quadratic source | share, uniform source | change")
    shares_q = {}
    for d, (t_pos, t_neg) in enumerate(((1, 2), (3, 4), (5, 6))):
        shp = (pw_q[t_pos] + pw_q[t_neg]) / sp
        shu = (pw_u[t_pos] + pw_u[t_neg]) / sp_u
        shares_q[d] = shp
        note = ""
        if d == 1:
            note = (f"  (y-pair share predicted ≈ 0.25, range [0.15, 0.40], from "
                    f"{REC['uni_y_share']:.6f}; in range: {0.15 <= shp <= 0.40}; not asserted)")
        say(f"  ±{AXNAME[d]}  {shp:.6f} | {shu:.6f} | {shp - shu:+.6f}{note}")
    clk.lap("face tags and powers")

    # ---- y centre line (x = z = 0), 129 points ---------------------------------------------
    yl = np.linspace(-0.995 * B, 0.995 * B, 129)
    pts = np.zeros((yl.size, 3))
    pts[:, 1] = yl
    t_q = np.array(_point_vals(sq.temperature, pts))
    t_uf = np.array(_point_vals(st.temperature, pts))
    t_se = np.real(g9._quad_ex()(pts[:, 0], pts[:, 1], pts[:, 2]))
    assert not np.any(np.isnan(t_q)) and not np.any(np.isnan(t_uf))
    j = int(np.argmax(t_q))
    say()
    say(f"[line] y centre line, 129 points: max theta (quadratic) {t_q[j]:.6f} K at |y|/B = "
        f"{abs(yl[j]) / B:.4f} (predicted in [0.75, 0.97]; in range: "
        f"{0.75 <= abs(yl[j]) / B <= 0.97}; not asserted)")
    k0 = int(np.argmin(np.abs(yl)))
    say(f"[line] theta at y = {yl[k0] / B:+.4f} B: quadratic {t_q[k0]:.6f}, series {t_se[k0]:.6f}, "
        f"uniform-source FEM {t_uf[k0]:.6f} K")
    say(f"[line] max |P2 - series| along the line {np.max(np.abs(t_q - t_se)):.4e} K; "
        f"max uniform-source theta {np.max(t_uf):.6f} K")
    clk.lap("centre line")

    # ---- centre-versus-side ratio at t* and steady ------------------------------------------
    r_t = tv[0] / tv[1]
    r_s = pv[0] / pv[1]
    say(f"[ratio] theta(0)/theta(0,0.9B,0): at t* = {T_STAR:.1f} s {r_t:.6f} (record "
        f"{REC['tv0'] / REC['tv09y']:.6f}); steady {r_s:.6f} (record "
        f"{REC['th0'] / REC['th09y']:.6f}); a probe at the centre reads ≈ {100 * r_t:.1f} % of the "
        "side-wall rise at t*; probe-placement reading only, no compliance claim")

    # ---- outputs: CSV, then XDMF, then adopt_host_ownership -----------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_y_centre_line.csv"
        with open(csv_path, "w") as fh:
            fh.write("y_m,y_over_B,theta_quad_K,theta_series_K,theta_uniform_K\n")
            for i in range(yl.size):
                fh.write(f"{yl[i]:.9e},{yl[i] / B:.6f},{t_q[i]:.9e},{t_se[i]:.9e},"
                         f"{t_uf[i]:.9e}\n")
    comm.Barrier()

    V1 = fem.functionspace(msh, ("Lagrange", 1))

    def p1(name, src, offset=0.0):
        f = fem.Function(V1, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - offset
        f.x.scatter_forward()
        return f

    f_q = fem.Function(V1, name="q_quad")  # P1 interpolant of the polynomial, display only
    f_q.interpolate(lambda x: Q * 3.0 * (x[1] ** 2 + x[2] ** 2) / (B_BOX**2 + C_BOX**2))
    f_q.x.array[:] = np.real(f_q.x.array)
    f_q.x.scatter_forward()
    f_p = p1("theta_quad_steady", sq.temperature, T_INF)
    f_se = p1("theta_quad_series_steady", ex_q)
    f_u = p1("theta_uniform_steady", st.temperature, T_INF)
    f_z = p1("theta_noz_series_steady", ex_z)
    f_t = p1("theta_tstar_quad", res.temperature, T_INF)
    f_ts = p1("theta_tstar_quad_series", ex_t)
    f_pd = fem.Function(V1, name="theta_quad_minus_series")
    f_pd.x.array[:] = np.real(f_p.x.array) - np.real(f_se.x.array)
    f_zd = fem.Function(V1, name="theta_quad_minus_noz_series")
    f_zd.x.array[:] = np.real(f_p.x.array) - np.real(f_z.x.array)
    f_pd.x.scatter_forward()
    f_zd.x.scatter_forward()
    n_own = V1.dofmap.index_map.size_local
    thin_local = float(np.max(np.abs(np.real(f_zd.x.array[:n_own])))) if n_own else 0.0
    thin = comm.allreduce(thin_local, op=MPI.MAX)
    say(f"[print] thin-axis term: max |theta_quad - theta_noz-series| over owned P1 dofs {thin:.6f} K "
        f"beside (c-z) {c_z:.4e} (record {REC['c_z']:.4e}, {GATE}:79)")
    fields = {}
    for f in (f_q, f_p, f_se, f_pd, f_u, f_z, f_zd, f_t, f_ts):
        fields[f.name] = f
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields, comm=comm,
        facet_tags=face_tags)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    clk.lap("outputs")
    if rank0:
        print(f"\n[output] {xdmf_path.name}: CG1 q_quad, theta_quad_steady, "
              "theta_quad_series_steady, theta_quad_minus_series, theta_uniform_steady, "
              "theta_noz_series_steady, theta_quad_minus_noz_series, "
              f"theta_tstar_quad, theta_tstar_quad_series (K, t* = {T_STAR:.1f} s), six face tags")
        print(f"[output] {BASENAME}_y_centre_line.csv ({yl.size} rows)")
        print(f"[time] total {time.perf_counter() - clk.t0:.1f} s")
        print("[scope] one analytic source q = 3Q(y^2 + z^2)/(B^2 + C^2) on the m = 4 rung at degree 2, "
              "-n 2; no EM solved, no coil, no h-rate, no F2182 compliance, no external check claimed")
        print("All assertions hold")


if __name__ == "__main__":
    main()
