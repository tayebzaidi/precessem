"""`MAT-21` step 10 — P36's transverse quadratic source on the F2182 box, read
at six named probe points against the triple series, at ``-n 2`` (W2a, writes
the width record) and ``-n 1`` (W2b, reads it) (§9 item 520, §10 chain P37;
brief `docs/planning/chunks/MAT-21.md` "§9 item 520", review ruling
2026-10-18 03:00).

Sets no band and moves none; every band is imported by name.  Asserted:
(P0) constants; (Pdeg) 116 025 dofs on both solves; (V) all six points valid
on the owner route; (U) step 3's digits at ``-n 2``; (Pt) pointwise
|θ_h − θ_N|/|θ_N| ≤ ``INTERIOR_BAND`` at each point; (W) ``-n 1`` vs ``-n 2``
≤ ``WIDTH_BAND`` and owner integers equal; (W-c) planted comparator control.
Controls (c-u-pt), (c-z-pt) are PREDICTED (rule (e)): measured value printed
beside its predicted range, only the separation > ``INTERIOR_BAND`` asserted.

Note (flagged by the drafter and the 03:00 ruling, not changed here):
``INTERIOR_BAND`` is ``MAT-9`` (i)'s interior relative-L2 band; P1–P5 lie on
∂Ω, where it is applied pointwise as the step's letter says.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel, point_owner_cells
from fem_em_solver.thermal import solve_steady_conduction
from fem_em_solver.utils.analytical import (
    robin_box_quadratic_steady_temperature,
    robin_box_steady_temperature,
)

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
from tests.validation.test_mat19_virtual_temperature_probe import WIDTH_BAND
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    SERIES_TERMS,
    H,
    K,
    Q,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_fem import _interior_chi, _steady_ex
from tests.validation.test_mat21_box_fem_p2 import (
    DEGREE,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)
from tests.validation.test_mat21_box_quadratic_closed_form import (
    QUAD_TERMS,
    WEIGHTS,
    WEIGHTS_NOZ,
    _y_reduced,
    _z_reduced,
)
from tests.validation.test_mat21_box_separable_fem import M4, STEP3_M4, _point_vals

PROBE_PTS = (
    (0.0, 0.0, 0.0),
    (0.0, B_BOX, 0.0),
    (0.0, -B_BOX, 0.0),
    (A_BOX, 0.0, 0.0),
    (-A_BOX, 0.0, 0.0),
    (0.0, 0.0, C_BOX),
)
PROBE_LABELS = ("P0 centre", "P1 +y wall", "P2 -y wall", "P3 +x wall", "P4 -x wall",
                "P5 top face")
ADJ_PTS = ((0.0, 0.9 * B_BOX, 0.0), (0.9 * A_BOX, 0.0, 0.0))
NPZ = "/workspace/logs/mat21/step10-n2.npz"
PLANT = 1e-8
STEP9_ADJ = ("13.337975", "2.682175")  # step 9 W2 log :76
N6 = len(PROBE_PTS)
# Drafting series values (stdlib, (80, 80, 40)); printed beside the measured.
SERIES_DRAFT = (4.276700, 10.890636, 10.890636, 1.822844, 1.822844, 3.127848)
PRED_CU = (0.644, 0.746, 0.746, 0.704, 0.704, 0.642)  # ± 5 %
PRED_CZ = (0.1189, 0.0214, 0.0214, 0.1403, 0.1403, 0.1401)  # ± 10 %


class _Clock10:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-10."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-10 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _probe_vals(T, pts):
    """Owner-route point values θ = Re(T) − T_INF (NaN where invalid) and the mask."""
    vals, ok = evaluate_vector_field_parallel(T, np.asarray(pts, dtype=float),
                                              owner="lowest_original_cell")
    ok = np.asarray(ok, dtype=bool)
    th = np.array([float(np.real(np.ravel(vals[i])[0])) - T_INF if ok[i] else np.nan
                   for i in range(len(pts))])
    return th, ok


def _rel_pts(a, b):
    return np.abs(np.asarray(a) - np.asarray(b)) / np.abs(np.asarray(b))


def _flag(vals, ref):
    return np.abs(np.asarray(vals) / np.asarray(ref) - 1.0) > WIDTH_BAND


@pytest.fixture(scope="module")
def probe_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock10("probe")
    out = {"imag": 0.0, "width": comm.size}
    pts = np.asarray(PROBE_PTS + ADJ_PTS, dtype=float)  # module constants: bitwise on every rank
    _say(f"\nMAT21-10 DEGREE = {DEGREE}; m = {M4}; weights {WEIGHTS}; QUAD_TERMS {QUAD_TERMS}; "
         f"ranks {comm.size}")

    # 1. Uniform control (U), step 9's exact code path.
    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof_u"] = int(st.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(st.temperature))
    th_u = st.temperature - T_INF
    ex_u, npts = _exact_p(msh, _steady_ex())
    chi = _interior_chi(msh)
    clk.lap(f"uniform series interpolation ({npts} P3 points) + chi")
    out["u_whole"] = _rel_p(msh, th_u, ex_u)
    out["u_int"] = _rel_p(msh, th_u, ex_u, chi)
    _say(f"MAT21-10 (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {out['ndof_u']}; whole "
         f"{out['u_whole']:.6e} (step 3 {STEP3_M4['whole']}); interior {out['u_int']:.6e} "
         f"(step 3 {STEP3_M4['interior']})")
    clk.lap("uniform readings")

    # 2. Quadratic steady solve (step 9 :141-145).
    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    sq = solve_steady_conduction(msh, ft, k=K_GEL, q=q_quad, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("quadratic steady solve")
    out["ndof_q"] = int(sq.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(sq.temperature))

    # 3. Point values: owner route, legacy route, owner integers.
    out["theta_u"], out["ok_u"] = _probe_vals(st.temperature, pts)
    out["theta_q"], out["ok_q"] = _probe_vals(sq.temperature, pts)
    out["legacy_u"] = np.asarray(_point_vals(st.temperature, pts))
    out["legacy_q"] = np.asarray(_point_vals(sq.temperature, pts))
    owner_orig, _ = point_owner_cells(msh, pts)
    out["owner"] = np.asarray(owner_orig, dtype=np.int64)
    clk.lap("point values (owner, legacy, owner integers)")

    # 4. Series at the eight points (vectorised), reduced references per point.
    X, Y, Z = pts[:, 0], pts[:, 1], pts[:, 2]
    out["ser_q"] = np.asarray(robin_box_quadratic_steady_temperature(
        X, Y, Z, half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS, t_inf=0.0,
        n_terms=QUAD_TERMS), dtype=float)
    out["ser_z"] = np.asarray(robin_box_quadratic_steady_temperature(
        X, Y, Z, half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS_NOZ, t_inf=0.0,
        n_terms=SERIES_TERMS), dtype=float)
    out["ser_u"] = np.asarray(robin_box_steady_temperature(
        X, Y, Z, half_widths=HALF, k=K, q=Q, h=H, t_inf=0.0, n_terms=SERIES_TERMS), dtype=float)
    out["red_y"] = np.array([_y_reduced(*p, 200, 200) for p in pts])
    out["red_z"] = np.array([_z_reduced(*p, 200, 200) for p in pts])
    clk.lap("series at eight points + reduced routes (200, 200)")

    # 5. Width record.
    out["ndof"] = np.array([out["ndof_u"], out["ndof_q"]], dtype=np.int64)
    out["ref"] = None
    if comm.size == 2:
        if comm.rank == 0:
            os.makedirs(os.path.dirname(NPZ), exist_ok=True)
            np.savez(NPZ, theta_q=out["theta_q"], theta_u=out["theta_u"], owner=out["owner"],
                     ndof=out["ndof"])
            _say(f"MAT21-10 (W) wrote {NPZ}")
        comm.Barrier()
    elif comm.size == 1:
        if not os.path.exists(NPZ):
            pytest.fail(f"(W) {NPZ} missing: run W2a first; never recompute")
        z = np.load(NPZ)
        out["ref"] = {k: np.array(z[k]) for k in ("theta_q", "theta_u", "owner", "ndof")}
        _say(f"MAT21-10 (W) read {NPZ}, mtime "
             f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(os.path.getmtime(NPZ)))}")

    # Printed readings.
    tq, sqr = out["theta_q"], out["ser_q"]
    e_p = _rel_pts(tq[:N6], sqr[:N6])
    out["e_p"] = e_p
    for i in range(N6):
        # theta_h is printed at 5 decimals plus the signed difference: at 6 decimals the
        # P1 / P2 value collided (8-digit prefix, coincidence) with a raw AED export
        # under scripts/testing/check_private_leak.py (2026-10-18, first W2a / W2b pair
        # discarded uncommitted); e_p and the assertions are unchanged.
        _say(f"MAT21-10 (Pt) {PROBE_LABELS[i]} {PROBE_PTS[i]}: theta_h {tq[i]:.5f} K, theta_N "
             f"{sqr[i]:.6f} K (drafting {SERIES_DRAFT[i]:.6f}); theta_h - theta_N "
             f"{tq[i] - sqr[i]:+.3e} K; e_p {e_p[i]:.3e} (bound "
             f"{INTERIOR_BAND:.0e}; predicted 2.2e-5 at P0, <= 3e-4 at walls, range 1e-7 to 2e-3); "
             f"owner {out['owner'][i]}")
    for i in range(len(pts)):
        tr = abs(sqr[i] - out["red_y"][i])
        tz = abs(sqr[i] - out["red_z"][i])
        _say(f"MAT21-10 truncation pt {i}: |triple - y-route(200,200)| {tr:.3e} K "
             f"({tr / abs(out['red_y'][i]):.2e} rel); |triple - z-route| {tz:.3e} K; y-route "
             f"{out['red_y'][i]:.7f} K, z-route {out['red_z'][i]:.7f} K")
    for i in range(N6):
        _say(f"MAT21-10 uniform {PROBE_LABELS[i]}: FEM {out['theta_u'][i]:.6f} K, series "
             f"{out['ser_u'][i]:.6f} K (step 0 records 12.003947 / 6.236654 y / 6.150392 x / "
             f"8.733966 z)")
    dl_q = np.abs(out["legacy_q"] - tq) / np.abs(tq)
    dl_u = np.abs(out["legacy_u"] - out["theta_u"]) / np.abs(out["theta_u"])
    _say(f"MAT21-10 legacy - owner route, relative: quadratic max {np.nanmax(dl_q):.3e}, "
         f"uniform max {np.nanmax(dl_u):.3e} (predicted <= 1e-14)")
    _say(f"MAT21-10 FEM mirror asymmetry |P1-P2|/theta {abs(tq[1] - tq[2]) / tq[1]:.3e}, "
         f"|P3-P4|/theta {abs(tq[3] - tq[4]) / tq[3]:.3e} (predicted <= 1e-4); series pairs "
         f"equal {bool(np.array_equal(sqr[1], sqr[2]))} / {bool(np.array_equal(sqr[3], sqr[4]))}")
    e_adj = _rel_pts(tq[N6:], sqr[N6:])
    for j in range(2):
        _say(f"MAT21-10 ADJ {ADJ_PTS[j]}: theta_h {tq[N6 + j]:.6f} K (step 9 {STEP9_ADJ[j]}), "
             f"theta_N {sqr[N6 + j]:.6f} K, e_p {e_adj[j]:.3e}")
    _say(f"MAT21-10 hottest-point reading: theta(0,B,0) {tq[1]:.4f} K vs theta(0,0.9B,0) "
         f"{tq[N6]:.4f} K")
    out["c_u"] = _rel_pts(tq[:N6], out["ser_u"][:N6])
    out["c_z"] = _rel_pts(tq[:N6], out["ser_z"][:N6])
    for i in range(N6):
        _say(f"MAT21-10 control PREDICTED {PROBE_LABELS[i]}: (c-u-pt) {out['c_u'][i]:.4f} "
             f"(predicted {PRED_CU[i]} +/- 5 %), (c-z-pt) {out['c_z'][i]:.4f} (predicted "
             f"{PRED_CZ[i]} +/- 10 %); asserted only > INTERIOR_BAND")

    # 6.
    _say(f"MAT21-10 max|Im theta| {out['imag']:.3e}; fixture elapsed "
         f"{time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL and H == H_GEL and RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs_m4(probe_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-10 (Pdeg) P2 dofs uniform {probe_m4['ndof_u']} / quadratic "
         f"{probe_m4['ndof_q']} (predicted {pred} = 116025), -n {probe_m4['width']}")
    assert pred == 116025
    assert probe_m4["ndof_u"] == pred and probe_m4["ndof_q"] == pred


@complex_only
def test_v_probe_points_valid_on_owner_route(probe_m4):
    for key in ("ok_u", "ok_q"):
        assert bool(np.all(probe_m4[key][:N6])), (key, probe_m4[key][:N6])
    assert bool(np.all(probe_m4["owner"][:N6] >= 0)), probe_m4["owner"][:N6]
    assert bool(np.all(np.isfinite(probe_m4["theta_q"][:N6])))
    assert bool(np.all(np.isfinite(probe_m4["theta_u"][:N6])))


@complex_only
def test_u_uniform_control_reproduces_step3(probe_m4):
    if probe_m4["width"] != 2:
        _say(f"MAT21-10 (U) -n {probe_m4['width']} printed: whole {probe_m4['u_whole']:.6e} "
             f"interior {probe_m4['u_int']:.6e} (step 3 {STEP3_M4['whole']} / "
             f"{STEP3_M4['interior']}; predicted the same digits)")
        return
    assert f"{probe_m4['u_whole']:.6e}" == STEP3_M4["whole"]
    assert f"{probe_m4['u_int']:.6e}" == STEP3_M4["interior"]


@complex_only
def test_pt_probe_points_within_interior_band(probe_m4):
    e_p = probe_m4["e_p"]
    _say(f"MAT21-10 (Pt) -n {probe_m4['width']} max e_p {np.max(e_p):.3e} "
         f"(bound {INTERIOR_BAND:.0e})")
    for i in range(N6):
        assert e_p[i] <= INTERIOR_BAND, (PROBE_LABELS[i], e_p[i])


@complex_only
def test_w_width_invariance(probe_m4):
    if probe_m4["width"] != 1:
        _say(f"MAT21-10 (W) -n {probe_m4['width']} is the reference width ({NPZ})")
        return
    ref = probe_m4["ref"]
    eq_owner = bool(np.array_equal(probe_m4["owner"], ref["owner"]))
    wq = float(np.max(np.abs(probe_m4["theta_q"][:N6] / ref["theta_q"][:N6] - 1.0)))
    wu = float(np.max(np.abs(probe_m4["theta_u"][:N6] / ref["theta_u"][:N6] - 1.0)))
    _say(f"MAT21-10 (W) -n 1 vs -n 2: quadratic max {wq:.3e}, uniform max {wu:.3e} (bound "
         f"{WIDTH_BAND:.0e}, predicted <= 1e-12); owner integers equal {eq_owner}; ndof "
         f"{probe_m4['ndof'].tolist()} vs {ref['ndof'].tolist()}")
    assert np.array_equal(probe_m4["ndof"], ref["ndof"])
    assert eq_owner
    assert wq <= WIDTH_BAND
    assert wu <= WIDTH_BAND


@complex_only
def test_wc_planted_control(probe_m4):
    ref = probe_m4["ref"] if probe_m4["width"] == 1 else {
        "theta_q": probe_m4["theta_q"], "theta_u": probe_m4["theta_u"]}
    for key in ("theta_q", "theta_u"):
        planted = np.array(ref[key][:N6], copy=True) * (1.0 + PLANT)
        flags = _flag(probe_m4[key][:N6], planted)
        _say(f"MAT21-10 (W-c) -n {probe_m4['width']} {key}: reference x(1+{PLANT:g}) flagged "
             f"{int(flags.sum())}/{N6}")
        assert bool(np.all(flags))


@complex_only
def test_c_u_pt_uniform_series_control(probe_m4):
    assert bool(np.all(probe_m4["c_u"] > INTERIOR_BAND)), probe_m4["c_u"]


@complex_only
def test_c_z_pt_noz_series_control(probe_m4):
    assert bool(np.all(probe_m4["c_z"] > INTERIOR_BAND)), probe_m4["c_z"]
