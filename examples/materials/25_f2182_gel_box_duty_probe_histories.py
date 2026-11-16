"""Example (`EX-131`, `mat:25`): the F2182-sized gel box under the duty-cycled transverse eddy-SAR
source, read as temperature histories at the six F2182 probe points.

`MAT-21` step 11 (W2) gates P36's source q = 3Q (y^2 + z^2) / (B^2 + C^2) switched by step 2's
schedule 11111000011110000111100001111000011110000 (ten switches at HALF_T = 282.15 s, duty 0.5)
on step 3's m = 4 rung at degree 2, CN N = 40 to t*, against a one-tensor ten-run closed form.
`MAT-21` step 10 gates the same source at six named probe points on the owner route.  Neither
gate shows a history.  Here the gate's three CN runs (default, chi = 1, duty) are re-run with
their own calls (copied; the pytest fixture cannot be called, and snapshots are added on the
chi = 1 and duty runs only), and theta(t_n) at the six probe points, every CN grid time t_n =
n * DT_BOX, n = 0 ... 40, is read through the owner route and printed beside the ten-run closed
form and the chi = 1 run with its series.

Everything -- every band, every helper -- is imported from the gate modules
(``tests.validation.test_mat21_box_quadratic_duty_fem`` as ``g11``,
``tests.validation.test_mat21_box_probe_points_fem`` as ``g10``,
``tests.validation.test_mat21_box_quadratic_duty_closed_form`` as ``w1``).  Records live only in
the example-local ``REC`` dict, each commented with its ``log:line``; they are printed beside the
readings and never asserted by equality.

Asserted (imported bands only): (Pdeg) dofs on the uniform, steady and three CN solves; (U) step
3's digits; (V) the six probes valid on the owner route; (Pt) steady pointwise error at each
probe; (S0) census and schedule string; (R36) the default run reproduces P36's digits; (D'-quad)
chi = 1 vs default; (T-qduty) duty CN interior vs the duty closed form; (E-qduty) balance on all
three CN runs.  Example-side identities, exact by construction: (H40) the n = 40 snapshot is the
run's own final field; (H=) the duty and chi = 1 snapshots coincide while both sources are on;
(D0-h) the one-tensor closed form at one switch equals the public transient at every grid time.
Controls: (c-step), (c-u-pt), (c-z-pt) backed, asserted > INTERIOR_BAND; (cH) predicted, only the
separation asserted.
Printed, not asserted: the per-probe histories, the history errors, the duty / chi = 1 ratios, the
count of local maxima per probe.
Not claimed: width invariance (step 10's (W) is provenance only), pointwise history accuracy as a
gate, other duty fractions, a periodic state, dose or the six-minute window, the box's EM SAR, a
coil, perfusion, a container wall, F2182 compliance or probe advice, any external check.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:25 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_25_f2182_gel_box_duty_probe_histories_history_combined.xdmf``,
``..._combined.xdmf`` and ``..._probe_histories.csv``.
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

from fem_em_solver.io.paraview_utils import (
    adopt_host_ownership,
    write_xdmf_time_series,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_mat21_box_probe_points_fem as g10  # noqa: E402
import tests.validation.test_mat21_box_quadratic_duty_closed_form as w1  # noqa: E402
import tests.validation.test_mat21_box_quadratic_duty_fem as g11  # noqa: E402

# Imported by name through the gate namespaces (never restated).
INTERIOR_BAND = g11.INTERIOR_BAND
D_BAND = g11.D_BAND
E_BAND = g11.E_BAND
ROBIN_TAG = g11.ROBIN_TAG
RHO_C = g11.RHO_C
K_GEL, H_GEL, N_MAIN, T_INF = g11.K_GEL, g11.H_GEL, g11.N_MAIN, g11.T_INF
Q, B_BOX, C_BOX = g11.Q, g11.B_BOX, g11.C_BOX
T_STAR, DT_BOX, HALF_T = g11.T_STAR, g11.DT_BOX, g11.HALF_T
DEGREE, M4, STEP3_M4 = g11.DEGREE, g11.M4, g11.STEP3_M4
SCHEDULE, R36 = g11.SCHEDULE, g11.R36
DUTY = g11.DUTY
QUAD_TERMS, WEIGHTS = g11.QUAD_TERMS, g11.WEIGHTS
_allsum, _max_imag = g11._allsum, g11._max_imag
_census, _chi_box = g11._census, g11._chi_box
_quad_trans_ex, _quad_duty, _quad_floor = g11._quad_trans_ex, g11._quad_duty, g11._quad_floor
_quad_duty_mean = g11._quad_duty_mean
_rung_p, _exact_p, _rel_p, _dofs_pred = g11._rung_p, g11._exact_p, g11._rel_p, g11._dofs_pred
_interior_chi, _steady_ex = g11._interior_chi, g11._steady_ex

PROBE_PTS, PROBE_LABELS = g10.PROBE_PTS, g10.PROBE_LABELS
_probe_vals, _rel_pts = g10._probe_vals, g10._rel_pts
WEIGHTS_NOZ, SERIES_TERMS, HALF, K, H = g10.WEIGHTS_NOZ, g10.SERIES_TERMS, g10.HALF, g10.K, g10.H
robin_box_quadratic_steady_temperature = g10.robin_box_quadratic_steady_temperature
robin_box_steady_temperature = g10.robin_box_steady_temperature
point_owner_cells = g10.point_owner_cells

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_25_f2182_gel_box_duty_probe_histories"

G11 = "20261018T124640Z_MAT-21-step11-W2-quadratic-duty-fem.log"
G10 = "20261018T111612Z_MAT-21-step10-W2a-points-n2.log"
W1 = "20261018T124327Z_MAT-21-step11-W1-closed-form.log"

# Records, printed beside readings only; never asserted by equality.
REC = {
    "u_whole": "2.489307e-05", "u_int": "1.669359e-05",             # G11 :70
    "r36_int": "2.694340e-04", "r36_whole": "4.310154e-04",         # G11 :74
    "dprime": 3.025e-15,                                            # G11 :76
    "t_qduty": 2.543310e-04, "t_whole": 3.319219e-04,               # G11 :82
    "e": 2.288e-13,                                                 # G11 :83 (round-off floor)
    "c_step": 0.50333,                                              # G11 :84
    "mean_duty": 0.598483, "mean_cf": 0.598518, "mean_chi1": 1.209940,  # G11 :85
    "drop_mean": 0.49464, "drop_centre": 0.54476, "drop_09y": 0.49470,  # G11 :86
    "floor": 1.2484e-06,                                            # G11 :88
    "e_p": (2.255e-05, 1.839e-05, 1.839e-05, 4.439e-05, 4.439e-05, 1.616e-04),  # G10 :76-81
    "theta_N": (4.2767, 10.890636, 10.890636, 1.822844, 1.822844, 3.127848),    # G10 :76-81
    "owner": (34782, 39774, 30103, 34938, 34636, 74718),            # G10 :76-81
    "c_u": (0.6437, 0.7462, 0.7462, 0.7036, 0.7036, 0.6418),        # G10 :101-106
    "c_z": (0.1188, 0.0214, 0.0214, 0.1403, 0.1403, 0.1403),        # G10 :101-106
    # W1 :76-81 closed-form duty/continuous ratio at t*: centre, (0,B,0) [= -y], (A,0,0), (0,0,C)
    "w1_ratio": (0.54466, 0.48560, 0.48560, 0.54159, 0.54159, 0.48572),
    "w1_mean_ratio": 0.49465,
}

N6 = len(PROBE_PTS)
NSTEP = N_MAIN
CHI_NAME = {"0": "off", "1": "on"}


class _Clock:
    """Per-phase wall clock, printed on rank 0 with flush, prefix EX-131."""

    def __init__(self):
        self.t = time.perf_counter()
        self.t0 = self.t

    def lap(self, phase):
        now = time.perf_counter()
        if MPI.COMM_WORLD.rank == 0:
            print(f"EX-131 clock {phase}: {now - self.t:.1f} s", flush=True)
        self.t = now


def _cn(msh, ft, q_quad, schedule, snaps=()):
    kw = {} if schedule is None else {"source_schedule": schedule}
    return solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_quad, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE, snapshot_steps=snaps, **kw)


def _n_maxima(h):
    """Local maxima of a history: strictly above the previous, not below the next."""
    return int(sum(1 for i in range(1, len(h) - 1) if h[i] > h[i - 1] and h[i] >= h[i + 1]))


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
    # g11's pin equalities (test_precondition_pins, :187-193)
    assert (g11.N_STEPS, g11.SPH_PIN, g11.N_SWITCH_PIN, g11.DUTY_PIN) == (
        N_MAIN, g11.STEPS_PER_HALF, g11.N_SWITCH, DUTY)
    assert T_STAR == 0.2 * g11.T_C
    assert T_STAR == g11.T_STAR_2
    assert T_INF == g11.T_INF_0
    assert HALF_T == g11.STEPS_PER_HALF * DT_BOX
    assert DT_BOX == T_STAR / N_MAIN
    assert DEGREE == 2 and g11.CMP_DEGREE == 3

    say("=" * 72)
    say("EX-131 -- the F2182-sized gel box, duty-cycled transverse eddy source: probe histories")
    say("=" * 72)
    say(f"[records] beside {G11}, {G10}, {W1}; DEGREE {DEGREE}, comparand P{g11.CMP_DEGREE}, "
        f"quadrature {g11.QDEG_P}, m = {M4}, weights {WEIGHTS}, QUAD_TERMS {QUAD_TERMS}, "
        f"ranks {comm.size}; DT_BOX {DT_BOX:.6f} s, HALF_T {HALF_T:.4f} s")
    say(f"[bands] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, D_BAND {D_BAND:.0e}, "
        f"E_BAND {E_BAND:.0e}")
    clk = _Clock()
    imag = 0.0
    ndof = {}
    pts = np.asarray(PROBE_PTS, dtype=float)  # module constants: bitwise on every rank
    X, Y, Z = pts[:, 0], pts[:, 1], pts[:, 2]

    # ---- (i) uniform control (U), the gate's code path --------------------------------------
    msh, ft, ns, n_cells, v_h, st = _rung_p(M4, clk, DEGREE)
    tdim = msh.topology.dim
    imap_c = msh.topology.index_map(tdim)
    n_cells_all = imap_c.size_local + imap_c.num_ghosts
    cells = np.arange(n_cells_all, dtype=np.int32)
    cell_tags = dmesh.meshtags(msh, tdim, cells, np.full(n_cells_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm)"},
        title="mat:25 -- F2182-sized gel box, duty-cycled transverse source: z = 0 slice",
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    clk.lap("setup figure (no-op unless FEM_EM_SETUP_FIGURES=1)")
    ndof["uniform"] = int(st.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(st.temperature))
    th_u = st.temperature - T_INF
    ex_u, npts = _exact_p(msh, _steady_ex())
    chi = _interior_chi(msh)
    clk.lap(f"uniform series interpolation ({npts} P3 points) + chi")
    u_whole = _rel_p(msh, th_u, ex_u)
    u_int = _rel_p(msh, th_u, ex_u, chi)
    say(f"[anchor] (U) m = {M4} {tuple(ns)}: cells {n_cells}; P2 dofs {ndof['uniform']}; whole "
        f"{u_whole:.6e} (step 3 {STEP3_M4['whole']}, record {REC['u_whole']}, {G11}:70); "
        f"interior {u_int:.6e} (step 3 {STEP3_M4['interior']}, record {REC['u_int']})")
    del ex_u

    # ---- (ii) quadratic steady solve at the six probe points (g10:149-153, :169-176) --------
    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    sq = solve_steady_conduction(msh, ft, k=K_GEL, q=q_quad, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("quadratic steady solve")
    ndof["steady"] = int(sq.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(sq.temperature))
    tq, ok_q = _probe_vals(sq.temperature, pts)
    owner_orig, _ = point_owner_cells(msh, pts)
    owner = np.asarray(owner_orig, dtype=np.int64)
    ser_q = np.asarray(robin_box_quadratic_steady_temperature(
        X, Y, Z, half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS, t_inf=0.0,
        n_terms=QUAD_TERMS), dtype=float)
    ser_z = np.asarray(robin_box_quadratic_steady_temperature(
        X, Y, Z, half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS_NOZ, t_inf=0.0,
        n_terms=SERIES_TERMS), dtype=float)
    ser_u = np.asarray(robin_box_steady_temperature(
        X, Y, Z, half_widths=HALF, k=K, q=Q, h=H, t_inf=0.0, n_terms=SERIES_TERMS), dtype=float)
    e_p = _rel_pts(tq, ser_q)
    c_u = _rel_pts(tq, ser_u)
    c_z = _rel_pts(tq, ser_z)
    clk.lap("steady point values, owner integers, series")
    say()
    for i in range(N6):
        say(f"[anchor] (Pt) {PROBE_LABELS[i]} {PROBE_PTS[i]}: steady theta_h {tq[i]:.5g} K, "
            f"theta_N {ser_q[i]:.5g} K (record {REC['theta_N'][i]:.5g}); e_p {e_p[i]:.3e} (bound "
            f"{INTERIOR_BAND:.0e}; record {REC['e_p'][i]:.3e}, {G10}:76-81); owner {owner[i]} "
            f"(record {REC['owner'][i]})")
    say(f"[anchor] (V) owner route valid at all six: {bool(np.all(ok_q))}; owner integers >= 0: "
        f"{bool(np.all(owner >= 0))}")
    for i in range(N6):
        say(f"[control] {PROBE_LABELS[i]}: (c-u-pt) {c_u[i]:.4f} (backed, record {REC['c_u'][i]}), "
            f"(c-z-pt) {c_z[i]:.4f} (backed, record {REC['c_z'][i]}; no factor bar) "
            f"[{G10}:101-106]; asserted only > {INTERIOR_BAND:.0e}")

    # ---- (iii) the three CN runs (g11:130-150), snapshots on chi1 and duty only -------------
    snaps = range(N_MAIN + 1)
    comm.Barrier()
    r_def = _cn(msh, ft, q_quad, None)
    clk.lap(f"CN N={N_MAIN} default")
    ndof["default"] = int(r_def.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(r_def.temperature))
    e_bal = {"default": float(r_def.max_balance_residual)}
    th_def = r_def.temperature - T_INF
    ex_t, _ = _exact_p(msh, _quad_trans_ex(T_STAR))
    clk.lap(f"continuous series interpolation ({QUAD_TERMS})")
    r36_int = _rel_p(msh, th_def, ex_t, chi)
    r36_whole = _rel_p(msh, th_def, ex_t)
    say()
    say(f"[anchor] (R36) default CN vs continuous quadratic series: interior {r36_int:.6e} "
        f"(P36 {R36['interior']}, record {REC['r36_int']}), whole {r36_whole:.6e} (P36 "
        f"{R36['whole']}) [{G11}:74]")

    r_one = _cn(msh, ft, q_quad, lambda t: 1.0, snaps)
    clk.lap(f"CN N={N_MAIN} chi1 (snapshots)")
    ndof["chi1"] = int(r_one.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(r_one.temperature))
    e_bal["chi1"] = float(r_one.max_balance_residual)
    th_one = r_one.temperature - T_INF
    dprime = _rel_p(msh, th_one, th_def)
    mean_one = float(np.real(_allsum(comm, th_one * ufl.dx(domain=msh)))) / v_h
    say(f"[anchor] (D'-quad) chi = 1 vs default at t*: whole rel L2 {dprime:.3e} (bound D_BAND "
        f"{D_BAND:.0e}; record {REC['dprime']:.3e}, {G11}:76)")
    del r_def

    r_d = _cn(msh, ft, q_quad, _chi_box, snaps)
    clk.lap(f"CN N={N_MAIN} duty (snapshots)")
    ndof["duty"] = int(r_d.temperature.function_space.dofmap.index_map.size_global)
    imag = max(imag, _max_imag(r_d.temperature))
    e_bal["duty"] = float(r_d.max_balance_residual)
    th_d = r_d.temperature - T_INF
    mean_d = float(np.real(_allsum(comm, th_d * ufl.dx(domain=msh)))) / v_h
    sched = "".join(str(int(x)) for x in r_d.source_multipliers)
    n_on, n_off = _census(r_d, DT_BOX, 2, DUTY * T_STAR, "quad duty")
    ex_d, npts_d = _exact_p(msh, lambda x, y, z: _quad_duty(x, y, z, T_STAR))
    clk.lap(f"duty series interpolation ({npts_d} P3 points, {QUAD_TERMS})")
    t_int = _rel_p(msh, th_d, ex_d, chi)
    t_whole = _rel_p(msh, th_d, ex_d)
    c_step = _rel_p(msh, th_d, ex_t, chi)
    clk.lap("duty readings")
    say(f"[anchor] (S0) schedule (2 half-steps, then 39 steps): {sched} (must equal step 2's "
        f"{SCHEDULE}; on {n_on}, off {n_off}) [{G11}:78, :81]")
    say(f"[anchor] (T-qduty) duty CN interior vs duty closed form at t* {t_int:.6e} (bound "
        f"{INTERIOR_BAND:.0e}; record {REC['t_qduty']:.6e}, {G11}:82); whole {t_whole:.6e} "
        f"(record {REC['t_whole']:.6e})")
    say("[anchor] (E-qduty) max balance residual: " + ", ".join(
        f"{k_} {v:.3e}" for k_, v in e_bal.items()) + f" (bound {E_BAND:.0e}; record "
        f"{REC['e']:.3e}, {G11}:83; round-off floor, band only)")
    say(f"[control] (c-step) duty CN vs continuous series interior {c_step:.4e} (backed, record "
        f"{REC['c_step']:.5f}, {G11}:84; asserted only > {INTERIOR_BAND:.0e}; ratio to band "
        f"{c_step / INTERIOR_BAND:.1f})")
    cf_mean = _quad_duty_mean(T_STAR)
    say(f"[print] volume means at t*: <theta_duty> {mean_d:.5g} K (record {REC['mean_duty']:.5g}, "
        f"closed form {cf_mean:.5g}, record {REC['mean_cf']:.5g}), <theta_chi1> {mean_one:.5g} K "
        f"(record {REC['mean_chi1']:.5g}); <theta_duty>/<theta_chi1> {mean_d / mean_one:.5f} "
        f"(record {REC['drop_mean']}, W1 closed form {REC['w1_mean_ratio']})")
    say(f"[print] comparand floor _quad_floor() {_quad_floor():.4e} (record {REC['floor']:.4e}, "
        f"{G11}:88)")

    # ---- (iv) histories at the six probes, every CN grid time ---------------------------------
    t_grid = np.arange(N_MAIN + 1) * DT_BOX
    hist_d = np.empty((N_MAIN + 1, N6))
    hist_c = np.empty((N_MAIN + 1, N6))
    ok_hist = True
    for n in range(N_MAIN + 1):
        v, ok = _probe_vals(r_d.snapshots[n], pts)
        hist_d[n] = v
        ok_hist = ok_hist and bool(np.all(ok))
        v, ok = _probe_vals(r_one.snapshots[n], pts)
        hist_c[n] = v
        ok_hist = ok_hist and bool(np.all(ok))
    clk.lap(f"{2 * (N_MAIN + 1)} probe calls on the snapshots")
    cf_d = np.zeros((N_MAIN + 1, N6))
    cf_c = np.zeros((N_MAIN + 1, N6))
    d0h = True
    for n in range(1, N_MAIN + 1):
        t = n * DT_BOX
        cf_d[n] = np.asarray(_quad_duty(X, Y, Z, t), dtype=float)
        cf_c[n] = np.asarray(_quad_duty(X, Y, Z, t, switch_times=(0.0,)), dtype=float)
        pub = np.asarray(w1._public(X, Y, Z, t), dtype=float)
        d0h = d0h and bool(np.array_equal(cf_c[n], pub))
    clk.lap(f"closed-form histories ({3 * N_MAIN} calls)")

    # example-side identities
    h40 = {}
    for lab, r in (("duty", r_d), ("chi1", r_one)):
        eq = bool(np.array_equal(r.snapshots[N_MAIN].x.array, r.temperature.x.array))
        h40[lab] = bool(comm.allreduce(eq, op=MPI.LAND))
    h_eq = {n: _rel_p(msh, r_d.snapshots[n] - T_INF, r_one.snapshots[n] - T_INF)
            for n in range(1, 5)}
    c_h = _rel_p(msh, r_d.snapshots[8] - T_INF, r_one.snapshots[8] - T_INF)
    say()
    say(f"[identity] (H40) n = 40 snapshot array_equal the run's temperature (AND over ranks): "
        f"duty {h40['duty']}, chi1 {h40['chi1']}")
    say("[identity] (H=) duty vs chi1 snapshot, first on-phase (both execute the scheduled form "
        "with s = 1): " + ", ".join(f"n={n} {v:.3e}" for n, v in h_eq.items()) +
        f" (bound D_BAND {D_BAND:.0e}; predicted 0)")
    say(f"[identity] (D0-h) _quad_duty(switch_times=(0.0,)) array_equal w1._public at the six "
        f"points for n = 1..{N_MAIN}: {d0h} (extends W1 (D0), {W1}:49-50, to the grid)")
    say(f"[control] (cH) duty vs chi1 snapshot at n = 8, end of the first off-phase: {c_h:.4f} "
        f"(PREDICTED 0.50, range [0.40, 0.60], ceiling ~1; in range: {0.40 <= c_h <= 0.60}; "
        f"asserted only > {INTERIOR_BAND:.0e}; ratio to band {c_h / INTERIOR_BAND:.1f})")
    say(f"[anchor] (V) histories valid on the owner route at all {2 * (N_MAIN + 1)} snapshots: "
        f"{ok_hist}")

    say()
    say("[history] duty FEM theta (K) at n = 4k, k = 0..10 (t = n * DT_BOX; chi at step n: "
        "SCHEDULE[n]); closed form in brackets")
    ks = list(range(0, N_MAIN + 1, 4))
    say("  n      t (s) chi | " + " | ".join(f"{lab[:2]}" for lab in PROBE_LABELS))
    for n in ks:
        chn = SCHEDULE[n] if n >= 1 else "-"
        row = " | ".join(f"{hist_d[n, i]:.5g} [{cf_d[n, i]:.5g}]" for i in range(N6))
        say(f"  {n:2d} {t_grid[n]:9.2f} {chn}   | {row}")
    say("[history] chi1 FEM theta (K) at n = 4k; series in brackets")
    for n in ks:
        row = " | ".join(f"{hist_c[n, i]:.5g} [{cf_c[n, i]:.5g}]" for i in range(N6))
        say(f"  {n:2d} {t_grid[n]:9.2f}   | {row}")

    err_d = np.max(np.abs(hist_d[1:] - cf_d[1:]), axis=0)
    err_c = np.max(np.abs(hist_c[1:] - cf_c[1:]), axis=0)
    say()
    for i in range(N6):
        say(f"[print] {PROBE_LABELS[i]}: max_n |theta_h - theta_cf| duty {err_d[i]:.3e} K, chi1 "
            f"{err_c[i]:.3e} K (predicted <= 1.5e-3 K, range [1e-6, 5e-3]; in range: "
            f"{1e-6 <= err_d[i] <= 5e-3 and 1e-6 <= err_c[i] <= 5e-3})")
    rel_t = np.abs(hist_d[N_MAIN] - cf_d[N_MAIN]) / np.abs(cf_d[N_MAIN])
    for i in range(N6):
        say(f"[print] {PROBE_LABELS[i]}: relative error at t* {rel_t[i]:.3e} (predicted <= 1e-3, "
            f"range [1e-6, 3e-3]; in range: {1e-6 <= rel_t[i] <= 3e-3})")
    ratio_fem = hist_d[N_MAIN] / hist_c[N_MAIN]
    ratio_cf = cf_d[N_MAIN] / cf_c[N_MAIN]
    for i in range(N6):
        rr = abs(ratio_fem[i] / REC["w1_ratio"][i] - 1.0)
        say(f"[print] {PROBE_LABELS[i]}: theta_duty/theta_chi1 at t* FEM {ratio_fem[i]:.5f}, "
            f"closed form {ratio_cf[i]:.5f}, W1 record {REC['w1_ratio'][i]} ({W1}:76-80); "
            f"relative gap to the record {rr:.2e} (predicted <= 1e-3; in range: {rr <= 1e-3})")
    for i in range(N6):
        say(f"[print] {PROBE_LABELS[i]}: theta_chi1(t*)/theta_steady {hist_c[N_MAIN, i] / tq[i]:.5f} "
            f"(steady {tq[i]:.5g} K)")
    n_max = [_n_maxima(hist_d[:, i]) for i in range(N6)]
    for i in range(N6):
        arg = [n for n in range(1, N_MAIN) if hist_d[n, i] > hist_d[n - 1, i]
               and hist_d[n, i] >= hist_d[n + 1, i]]
        say(f"[print] {PROBE_LABELS[i]}: local maxima of the duty history {n_max[i]} at n = {arg} "
            "(drafter physics, no logged backing: P1/P2/P5 sawtooth with 5 maxima at "
            "n = 4, 12, 20, 28, 36; P0/P3/P4 monotone; printed only)")
    say(f"[print] max|Im theta| {imag:.3e}")

    pred = _dofs_pred(M4)
    say(f"[anchor] (Pdeg) P2 dofs {ndof} (predicted {pred})")

    # ---- assertions (imported bands; example-side identities; controls) ---------------------
    assert comm.size == 2 and DEGREE == 2 and M4 == 4
    assert pred == 116025, pred
    assert len(ndof) == 5 and all(v == pred for v in ndof.values()), f"(Pdeg) {ndof} vs {pred}"
    assert f"{u_whole:.6e}" == STEP3_M4["whole"], f"(U) whole {u_whole:.6e}"
    assert f"{u_int:.6e}" == STEP3_M4["interior"], f"(U) interior {u_int:.6e}"
    assert bool(np.all(ok_q)) and bool(np.all(owner >= 0)) and ok_hist, "(V)"
    assert bool(np.all(np.isfinite(tq))), "(V) finite"
    for i in range(N6):
        assert e_p[i] <= INTERIOR_BAND, f"(Pt) {PROBE_LABELS[i]} {e_p[i]:.3e}"
    assert sched == SCHEDULE, f"(S0) {sched}"
    assert f"{r36_int:.6e}" == R36["interior"], f"(R36) interior {r36_int:.6e}"
    assert f"{r36_whole:.6e}" == R36["whole"], f"(R36) whole {r36_whole:.6e}"
    assert dprime <= D_BAND, f"(D'-quad) {dprime:.3e}"
    assert t_int <= INTERIOR_BAND, f"(T-qduty) {t_int:.4e}"
    assert len(e_bal) == 3
    for name, r_ in e_bal.items():
        assert r_ <= E_BAND, f"(E-qduty) {name}: {r_:.3e}"
    assert h40["duty"] and h40["chi1"], f"(H40) {h40}"
    for n, v in h_eq.items():
        assert v <= D_BAND, f"(H=) n={n}: {v:.3e}"
    assert d0h, "(D0-h)"
    assert c_step > INTERIOR_BAND, f"(c-step) did not miss: {c_step:.4e}"
    assert bool(np.all(c_u > INTERIOR_BAND)), f"(c-u-pt) did not miss: {c_u}"
    assert bool(np.all(c_z > INTERIOR_BAND)), f"(c-z-pt) did not miss: {c_z}"
    assert c_h > INTERIOR_BAND, f"(cH) did not separate: {c_h:.4e}"

    # ---- outputs: frames, steady/state fields, CSV; adopt_host_ownership last ----------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    V1 = fem.functionspace(msh, ("Lagrange", 1))

    def p1(name, src, offset=0.0):
        f = fem.Function(V1, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - offset
        f.x.scatter_forward()
        return f

    frames = []
    for n in range(N_MAIN + 1):
        frames.append((float(n * DT_BOX), {
            "theta_duty": p1("theta_duty", r_d.snapshots[n], T_INF),
            "theta_chi1": p1("theta_chi1", r_one.snapshots[n], T_INF),
        }))
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_history_combined", msh, cell_tags, frames, comm=comm)
    clk.lap("history frames")

    f_se = p1("theta_duty_series_tstar", ex_d)
    f_dt = p1("theta_duty_tstar", r_d.temperature, T_INF)
    f_st = p1("theta_quad_steady", sq.temperature, T_INF)
    f_df = fem.Function(V1, name="theta_duty_minus_series_tstar")
    f_df.x.array[:] = np.real(f_dt.x.array) - np.real(f_se.x.array)
    f_df.x.scatter_forward()
    fields = {f.name: f for f in (f_st, f_dt, f_se, f_df)}
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields, comm=comm)

    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_probe_histories.csv"
        with open(csv_path, "w") as fh:
            fh.write("n,t_s,probe,x_m,y_m,z_m,theta_duty_fem_K,theta_duty_closed_K,"
                     "theta_chi1_fem_K,theta_chi1_series_K,chi\n")
            for n in range(N_MAIN + 1):
                chn = SCHEDULE[n] if n >= 1 else "-"
                for i in range(N6):
                    p = PROBE_PTS[i]
                    fh.write(f"{n},{t_grid[n]:.6f},{PROBE_LABELS[i].replace(' ', '_')},"
                             f"{p[0]:.6f},{p[1]:.6f},{p[2]:.6f},{hist_d[n, i]:.9e},"
                             f"{cf_d[n, i]:.9e},{hist_c[n, i]:.9e},{cf_c[n, i]:.9e},{chn}\n")
    comm.Barrier()
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    clk.lap("outputs")
    if rank0:
        print(f"\n[output] {series_path.name}: 41 frames (t = n * {DT_BOX:.4f} s) of P1 theta_duty "
              "and theta_chi1 (K), interpolated from the snapshots")
        print(f"[output] {xdmf_path.name}: CG1 theta_quad_steady, theta_duty_tstar, "
              f"theta_duty_series_tstar, theta_duty_minus_series_tstar (K, t* = {T_STAR:.1f} s)")
        print(f"[output] {BASENAME}_probe_histories.csv ({(N_MAIN + 1) * N6} rows)")
        print(f"[time] total {time.perf_counter() - clk.t0:.1f} s")
        print("[scope] one analytic source q = 3Q(y^2 + z^2)/(B^2 + C^2) under one duty schedule, "
              "m = 4, degree 2, -n 2, CN to t*, six points; no EM solved, no coil, no width "
              "invariance, no dose, no F2182 compliance, no external check claimed")
        print("All assertions hold")


if __name__ == "__main__":
    main()
