"""Example (`EX-69`): perfusion carries half the heat — the Pennes sink on the
Robin ball, steady and rising.

No other example has a perfusion term (`mat:4` and `mat:5` exclude it by name).
This one takes the 10 cm Robin-cooled ball of `MAT-11` / `MAT-13` (k = 0.6,
q = 2000 W/m^3, h = 10 W/m^2/K, arterial = ambient = 20 C) and switches the
sink w (theta - theta_a) on at w = 240 W/m^3/K (m R = 2):

* the **steady** field, perfused and unperfused — `MAT-11`; perfusion takes
  4.25 of the 8.35 W;
* the **rise** to t* = 0.2 R^2 rho c / k, Crank-Nicolson with two
  backward-Euler start-up half-steps, N = 40 — `MAT-13`; at t* the power splits
  storage / convection / perfusion = 21 / 39 / 40 %.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The fixture,
helpers and bands are imported from
``tests/validation/test_mat11_perfusion_sphere.py`` and
``tests/validation/test_mat13_transient_perfusion_sphere.py``; `MAT-13`'s
by-value copies of ``W``, ``INTERIOR_BAND`` and ``BALANCE_BAND`` are asserted
equal to `MAT-11`'s, where they are defined. No gate edit.

Asserted, each as its gate writes it: steady (i) interior relative L2 on theta
<= ``INTERIOR_BAND`` and (iii) balance <= ``BALANCE_BAND``; transient (P_w)
<= ``INTERIOR_BAND`` and (E) <= ``BALANCE_BAND``. Negative controls, asserted
only as "misses ``INTERIOR_BAND``": steady w = 0, transient w = 0.

Errors are on theta = T - T_inf, never on T (both gates run at T_inf = 20).

Run it::

    ./run_examples.sh -e mat:6 -n 2 -t 180

Output (on the ball mesh): ``examples/materials/paraview_output/
materials_06_perfusion_sink_robin_ball_steady_combined.xdmf`` (P1
``theta_perfused``, ``theta_unperfused``), ``..._transient_combined.xdmf`` (P1
``theta`` at steps 0, 4, ..., 40, w = 240) and ``..._rise_curve.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_time_series,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_steady_conduction  # noqa: E402

from tests.validation import test_mat11_perfusion_sphere as mat11  # noqa: E402
from tests.validation import test_mat13_transient_perfusion_sphere as mat13  # noqa: E402
from tests.validation.test_mat11_perfusion_sphere import (  # noqa: E402
    BALANCE_BAND,
    INTERIOR_BAND,
    W,
    _balance,
    _centre_surface,
    _errors,
)
from tests.validation.test_mat13_transient_perfusion_sphere import (  # noqa: E402
    H,
    K,
    Q,
    R,
    RHO_C,
    ROBIN_TAG,
    T_END,
    T_INF,
    _allsum,
    _ball,
    _exact_theta,
    _interior_error,
    _series_0,
    _series_w,
    _transient,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_06_perfusion_sink_robin_ball"
N_STEPS = 40
FRAMES = tuple(range(0, N_STEPS + 1, 4))

#: Gate records, complex build — `20261001T004409Z_MAT-11-cx-W2.log:119,130`
#: and `20261001T004517Z_MAT-13-cx-W3.log:99-102`. Printed beside this run, a
#: difference reported and never asserted.
RECORD = {
    "steady_i": 1.8285e-03, "steady_ctl": 1.1055e00,
    "conv_w": 4.095893, "perf_w": 4.251127, "src_w": 8.347021,
    "centre": 5.5505, "surface": 3.2660,
    "P_w": 1.6426e-03, "ctl_t": 3.8871e-01,
    "perf_t": 3.3197, "src_t": 8.3470,
    "split": (0.2134, 0.3889, 0.3977), "w0": 1.7085e-03,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _rec(key, value, fmt=".4e"):
    return f"{value:{fmt}} (record {RECORD[key]:{fmt}}, diff {value - RECORD[key]:+.1e})"


def _theta_of(T_h, name):
    """theta = T - T_inf as a new function in T_h's space."""
    f = fem.Function(T_h.function_space, name=name)
    f.x.array[:] = T_h.x.array - T_INF
    f.x.scatter_forward()
    return f


def _mean(comm, msh, f, vol):
    return float(np.real(_allsum(comm, f * ufl.dx(domain=msh)))) / vol


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # MAT-13 restates these by value; MAT-11 defines them. Structural check.
    assert mat13.W == W and mat13.INTERIOR_BAND == INTERIOR_BAND
    assert mat13.BALANCE_BAND == BALANCE_BAND
    assert (mat11.R, mat11.K, mat11.Q, mat11.H, mat11.T_INF) == (R, K, Q, H, T_INF)

    _say(comm, "=" * 78)
    _say(comm, "EX-69 -- perfusion carries half the heat: the Pennes sink on the Robin ball")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] bands imported from MAT-11: interior {INTERIOR_BAND:.0e}, "
               f"balance {BALANCE_BAND:.0e}; w = {W:g} W/m^3/K (m R = 2), "
               f"t* = {T_END:.4f} s, N = {N_STEPS}")

    t = time.perf_counter()
    msh, cell_tags, ft = _ball()
    n_cells = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] MAT-13's _ball(): {n_cells} cells")

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel ball (Pennes sink w = 240)"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:6 -- 10 cm Robin ball with perfusion (MAT-11 / MAT-13 fixture), h = 10",
        comm=comm,
    )

    # ---- steady (MAT-11) -----------------------------------------------------
    sol_w = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
        perfusion_w_per_m3_k=W,
    )
    sol_0 = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
    )
    e_in, e_all, vol, _n = _errors(msh, sol_w.temperature, W)
    bal = _balance(sol_w)
    centre, surf = _centre_surface(msh, sol_w.temperature)
    th_w = _theta_of(sol_w.temperature, "theta_perfused")
    th_0 = _theta_of(sol_0.temperature, "theta_unperfused")
    steady_mean = _mean(comm, msh, th_w, vol)
    steady_mean_0 = _mean(comm, msh, th_0, vol)
    e_ctl_s, _, _, _ = _errors(msh, sol_0.temperature, W)
    _say(comm, f"\n[steady] (i) interior rel L2 on theta, w = {W:g}: {_rec('steady_i', e_in)}  "
               f"ASSERTED <= {INTERIOR_BAND:.0e}")
    _say(comm, f"[steady] (iii) balance |conv + perf - src|/src = {bal:.3e}  "
               f"ASSERTED <= {BALANCE_BAND:.0e}")
    _say(comm, f"[steady] convected {sol_w.convected_power:.6f} W "
               f"(record {RECORD['conv_w']:.6f}), perfused {sol_w.perfused_power:.6f} W "
               f"(record {RECORD['perf_w']:.6f}), source {sol_w.source_power:.6f} W "
               f"(record {RECORD['src_w']:.6f}); perfusion carries "
               f"{sol_w.perfused_power / sol_w.source_power:.4f} of the source")
    _say(comm, f"[steady] theta_h centre {centre:.4f} K (record {RECORD['centre']:.4f}), "
               f"surface mean {surf:.4f} K (record {RECORD['surface']:.4f}); "
               f"volume mean {steady_mean:.4f} K perfused vs {steady_mean_0:.4f} K unperfused")
    _say(comm, f"[control] steady w = 0 against the perfused closed form: "
               f"{_rec('steady_ctl', e_ctl_s)}  ASSERTED > {INTERIOR_BAND:.0e}")

    # ---- transient (MAT-13), T_inf = T(0) = 20 ----------------------------------
    s = _transient(msh, ft, N_STEPS, snapshot_steps=FRAMES + (N_STEPS - 1,))
    s0 = _transient(msh, ft, N_STEPS, w=0.0, snapshot_steps=FRAMES)
    theta_ex = _exact_theta(msh, _series_w)
    e_pw = _interior_error(msh, s.temperature, theta_ex)
    e_ctl_t = _interior_error(msh, s0.temperature, theta_ex)
    e_w0 = _interior_error(msh, s0.temperature, _exact_theta(msh, _series_0))
    _say(comm, f"\n[transient] CN + Rannacher, N = {N_STEPS}, t* = {T_END:.4f} s")
    _say(comm, f"[transient] (P_w) interior rel L2 at t*: {_rec('P_w', e_pw)}  "
               f"ASSERTED <= {INTERIOR_BAND:.0e}")
    _say(comm, f"[transient] (E) max per-step balance {s.max_balance_residual:.3e}  "
               f"ASSERTED <= {BALANCE_BAND:.0e}")
    _say(comm, f"[transient] perfused power {s.perfused_power:.4f} W of {s.source_power:.4f} W "
               f"(record {RECORD['perf_t']:.4f} of {RECORD['src_t']:.4f}); "
               f"steady perfused {sol_w.perfused_power:.4f} W")

    # last-step power split, recomputed from snapshots 39 / 40 (as the gate)
    T39, T40 = s.snapshots[N_STEPS - 1], s.snapshots[N_STEPS]
    dt = T_END / N_STEPS
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    stor = float(np.real(_allsum(comm, RHO_C * (T40 - T39) / dt * ufl.dx)))
    conv = float(np.real(_allsum(comm, H * (0.5 * T40 + 0.5 * T39 - T_INF) * ds(ROBIN_TAG))))
    perf = float(np.real(_allsum(comm, W * (0.5 * T40 + 0.5 * T39 - T_INF) * ufl.dx)))
    src = s.source_power
    bal_E = abs(stor + conv + perf - src) / src
    _say(comm, f"[transient] (E) recomputed last-step balance {bal_E:.3e}  "
               f"ASSERTED <= {BALANCE_BAND:.0e}; split storage / convected / perfused = "
               f"{stor / src:.4f} / {conv / src:.4f} / {perf / src:.4f} "
               f"(record {RECORD['split'][0]:.4f} / {RECORD['split'][1]:.4f} / "
               f"{RECORD['split'][2]:.4f}; printed)")
    _say(comm, f"[control] transient w = 0 against the perfused series: "
               f"{_rec('ctl_t', e_ctl_t)}  ASSERTED > {INTERIOR_BAND:.0e}; "
               f"(w -> 0) against MAT-9's series {_rec('w0', e_w0)} (printed)")

    # ---- the rise curve (printed and written, never gated) -----------------------
    rows = []
    _say(comm, f"\n[rise] step   t [s]    mean theta_h [K] w=240   w=0        "
               f"steady mean w=240 {steady_mean:.4f} K  (printed, not gated)")
    for step in FRAMES:
        m_w = _mean(comm, msh, _theta_of(s.snapshots[step], "theta"), vol)
        m_0 = _mean(comm, msh, _theta_of(s0.snapshots[step], "theta"), vol)
        rows.append((step, step * dt, m_w, m_0))
        _say(comm, f"[rise] {step:4d}  {step * dt:8.3f}   {m_w:.6e}       {m_0:.6e}")

    # ---- assertions ---------------------------------------------------------------
    assert e_in <= INTERIOR_BAND, f"(i) {e_in:.4e}"
    assert bal <= BALANCE_BAND, f"(iii) {bal:.3e}"
    assert e_pw <= INTERIOR_BAND, f"(P_w) {e_pw:.4e}"
    assert s.max_balance_residual <= BALANCE_BAND, f"(E) {s.max_balance_residual:.3e}"
    assert bal_E <= BALANCE_BAND, f"(E) recomputed {bal_E:.3e}"
    assert e_ctl_s > INTERIOR_BAND, f"steady control {e_ctl_s:.4e} did not miss"
    assert e_ctl_t > INTERIOR_BAND, f"transient control {e_ctl_t:.4e} did not miss"

    # ---- outputs (on the ball mesh; CG1 exports are interpolants, nothing asserted) --
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    steady_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_steady_combined", msh, None,
        {"theta_perfused": th_w, "theta_unperfused": th_0}, comm=comm,
    )
    frames = [(step * dt, {"theta": _theta_of(s.snapshots[step], "theta")}) for step in FRAMES]
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_transient_combined", msh, None, frames, comm=comm,
    )
    table = OUTPUT_DIR / f"{BASENAME}_rise_curve.csv"
    if comm.rank == 0:
        with table.open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["step", "t_s", "mean_theta_w240_K", "mean_theta_w0_K",
                         "steady_mean_w240_K", "steady_mean_w0_K"])
            for step, tt, m_w, m_0 in rows:
                wr.writerow([step, f"{tt:.6f}", f"{m_w:.9e}", f"{m_0:.9e}",
                             f"{steady_mean:.9e}", f"{steady_mean_0:.9e}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {steady_path} -- P1 theta_perfused, theta_unperfused (theta = T - T_inf)")
    _say(comm, f"[output] {series_path} -- P1 theta (w = 240) at steps {FRAMES[0]}, {FRAMES[1]}, ..., {FRAMES[-1]}")
    _say(comm, f"[output] {table}")

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
