"""Example (`EX-110`): the six-minute window on an FEM heating history, and on a field.

`mat:17` windows closed forms; `mat:15` shows this coil-driven duty history
unwindowed.  This example joins them on `MAT-17` step 3's fixture (`MAT-15`
step 1's ``_run_heated``: the F-small coil-driven gel phantom at 10 MHz, drive
k = 0, G = 1 W / P_1, P = t*/5, D = 0.5): `MAT-20` step 2's FEM clause, re-run,
then the same operator applied to every P1 vertex of the phantom.

Nothing is restated.  ``_fem_history(keep_duty_fields=True)`` and
``_identities`` are `MAT-20` step 2's fixture body and identity computations
lifted to module level (rule (a) lift, EX-110; additive, the gate re-run green
in-slot); every band is imported from the gate modules.

Asserted in-script, each printed beside ``20261011T172810Z_MAT-20-step2-fem-W1.log``:
(H0) Duhamel on the per-step volume mean <= ``L_BAND``; (FEM) and (FEM-2d) the
six-minute window of the duty series equals the signed sum of the ten switch-on
runs' windows <= ``D_BAND``; control (Wc-FEM) shifted by P/2 > ``D_BAND``.
Example-side: (Cm) the volume mean of the per-vertex six-minute field equals the
windowed series <= ``D_BAND`` (both maps are linear); control (Cm-c) the window
one sample short misses.  Printed, never asserted: (Dc), the records, the power
against temperature windows, the early-peaking vertices, and Part B (the tuned
F-human's 10 g ball-mean CEM43 and six-minute peak from ``logs/mat15_s5/``).

Scope: operator identities on one FEM history; 1 W normalisation (no absolute
temperature), degree 1, no perfusion.  Not a limit, compliance, tissue, a periodic
state or a CEM43 gate.

Run through the harness (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:18 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_18_six_minute_window_coil_heating_series_combined.xdmf``,
``..._combined.xdmf`` and ``..._window_history.csv``.
"""

from __future__ import annotations

import json
import math
import sys
import types
from pathlib import Path

import numpy as np
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import (
    adopt_host_ownership,
    write_xdmf_time_series,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.thermal.dose import cem43, running_average

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_mat15_coil_driven_heating_phantom as m15  # noqa: E402
import tests.validation.test_mat17_duty_cycle_coil_driven as m17s3  # noqa: E402
import tests.validation.test_mat20_fem_history as m20f  # noqa: E402
import tests.validation.test_mat20_running_average as s2  # noqa: E402
import tests.validation.test_mat20_tuned_heating_dose as m20d  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_18_six_minute_window_coil_heating"
W1 = "20261011T172810Z_MAT-20-step2-fem-W1.log"
W2 = "20261011T183717Z_MAT-20-step3-W2.log"

# Record digits, example-local only (each beside its W1 / W2 line).
REC = {
    "h0": 5.689e-14,        # W1:1968
    "fem": 3.434e-14,       # W1:1969
    "fem2d": 3.444e-14,     # W1:1971
    "wc": 1.326e-01,        # W1:1972
    "dc": 0.043296,         # W1:1973
    "duty_32_36_40": (0.466250330, 0.584438780, 0.568835514),  # W1:1967
    "a_mx": 5.270483664e-01,   # W1:1974
    "h0_mx": 1.004390808e+00,  # W1:1974
    "t_max": 1254.0,           # W1:1974
    "h_tstar": (1.153313650151e+00, 1.050769107419e+00, 9.456669333318e-01,
                8.379629153190e-01, 7.275636005745e-01, 6.143603575039e-01,
                4.982231179698e-01, 3.789902177675e-01, 2.564498589888e-01,
                1.302990494639e-01),  # W1:1957-1966
    "b_cem43_s": 1.846214193739186e-11,   # W2:52 (cem43)
    "b_quad_s": 1.846214193739185e-11,    # W2:52 (quad)
    "b_peak_k": 4.277005459539476e-02,    # W2:57
    "b_dose_min": 3.077024e-13,           # W2:62
}


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    structural = {
        "m20f.D_BAND is m17s3.D_BAND": m20f.D_BAND is m17s3.D_BAND,
        "m20f.L_BAND is m17s3.L_BAND": m20f.L_BAND is m17s3.L_BAND,
        "m20f._run_heated is m15._run_heated": m20f._run_heated is m15._run_heated,
        "m20f.T_STAR is m15.T_STAR": m20f.T_STAR is m15.T_STAR,
        "m20f.W == s2.W_B2 == 360.0": m20f.W == s2.W_B2 == 360.0,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    D_BAND, L_BAND, DT, HALF = m17s3.D_BAND, m17s3.L_BAND, m17s3.DT, m17s3.HALF
    W, T_INF, T_STAR = m20f.W, m15.T_INF, m15.T_STAR
    if rank0:
        print("=" * 72)
        print("EX-110 -- the six-minute window on an FEM heating history (MAT-20 FEM clause)")
        print("=" * 72)
        print("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
        print(f"[records] Part A beside {W1}; Part B beside {W2}", flush=True)

    # ---- the gate's own fixture, with its fields ------------------------------
    r = m20f._fem_history(keep_duty_fields=True)
    base, sub, V1 = r["base"], r["sub"], r["V1"]
    work, f_mean, v_sub, q0 = r["work"], r["f_mean"], r["v_sub"], r["q0"]
    t, duty, hs = r["t"], r["duty"], r["H"]

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        base["msh"],
        base["cell_tags"],
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            CONDUCTOR_CELL_TAG: "conductor",
            2: "air",
            m15.PHANTOM_CELL_TAG: "phantom (gel)",
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(2,),
        translucent_tags=(m15.PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="mat:18 -- the six-minute window on the duty-cycled coil-heated gel phantom "
              "(MAT-20 FEM clause)",
        comm=comm,
    )

    ids = m20f._identities(t, duty, hs)
    ok = ids["ok"]
    a, scale = ids["a"], ids["scale"]

    # ---- (Cm) and (Cm-c): the window of the field -----------------------------
    F = np.stack(r["duty_fields"])  # (41, n_local + ghosts), theta = T - T_INF, real
    n_own = V1.dofmap.index_map.size_local

    def _volume_means_of_windowed_field(window):
        means, mx, tmx = running_average(t, F, window=window)
        fns = {}
        for n in m20f.STEPS:
            f = fem.Function(V1, name=f"theta_win_{n}")
            f.x.array[:] = T_INF + np.where(np.isnan(means[n]), 0.0, means[n])
            f.x.scatter_forward()
            fns[n] = f
        series = m20f._series(types.SimpleNamespace(snapshots=fns), work, f_mean, comm, v_sub)
        return means, mx, tmx, series

    means_F, mx_F, tmx_F, vm_F = _volume_means_of_windowed_field(W)
    cm = float(np.max(np.abs(vm_F[ok] - a[ok])) / scale)
    _, _, _, vm_short = _volume_means_of_windowed_field(W - DT)
    cmc = float(np.max(np.abs(vm_short[ok] - a[ok])) / scale)

    # ---- printed arithmetic: the source's own six-minute on-fraction ----------
    period = 2.0 * HALF

    def _on_time(tt):  # cumulative on-time of chi, exact: chi = 1 on [jP, jP + P/2)
        tt = max(tt, 0.0) + 1e-9
        return math.floor(tt / period) * HALF + min(tt % period, HALF)

    frac = np.full(t.size, np.nan)
    for n in np.nonzero(ok)[0]:
        frac[n] = (_on_time(t[n]) - _on_time(t[n] - W)) / W
    n_frac = int(np.nanargmax(frac))
    ratio = ids["a_mx"] / ids["h0_mx"]

    # ---- per-vertex statistics over owned dofs ---------------------------------
    own_mx = mx_F[:n_own]
    own_tmx = tmx_F[:n_own]
    owned = comm.allgather((int(n_own), int(F.shape[1])))
    th_tstar_max = comm.allreduce(float(np.max(F[-1][:n_own], initial=-np.inf)), op=MPI.MAX)
    mx_max = comm.allreduce(float(np.max(own_mx, initial=-np.inf)), op=MPI.MAX)
    n_tot = comm.allreduce(int(n_own), op=MPI.SUM)
    t36 = 36 * DT
    n_tstar = comm.allreduce(int(np.sum(np.isclose(own_tmx, T_STAR, atol=1e-6))), op=MPI.SUM)
    n_36 = comm.allreduce(int(np.sum(np.isclose(own_tmx, t36, atol=1e-6))), op=MPI.SUM)
    n_early = n_tot - n_tstar

    # ---- Part B (printed only) ---------------------------------------------------
    partb = None
    if m20d.SERIES_JSON.exists():
        with open(m20d.SERIES_JSON) as fh:
            d = json.load(fh)
        tb = np.asarray(d["times_s"], dtype=float)
        thb = np.asarray(d["theta_ball10g_k"], dtype=float)
        Tb = d["t_inf_c"] + thb
        dose = cem43(tb, Tb)
        quad = m20d._quad_dose(tb, Tb, math.log(m20d.R_BELOW))
        _, b_mx, b_tmx = running_average(tb, thb, window=W)
        partb = (dose, quad, b_mx, b_tmx, dose / 60.0)

    # ---- prints (rank 0) ----------------------------------------------------------
    if rank0:
        print(f"\n[anchor] (H0) {ids['h0']:.3e} <= {L_BAND:.0e}  (W1:1968 {REC['h0']:.3e})")
        print(f"[anchor] (FEM) W = {W} s over {int(ok.sum())} windows {ids['fem']:.3e} <= "
              f"{D_BAND:.0e}  (W1:1969 {REC['fem']:.3e})")
        print(f"[anchor] (FEM-2d) {ids['fem2d']:.3e} <= {D_BAND:.0e}; 2-D duty column == 1-D: "
              f"{ids['fem2d_eq']}  (W1:1971 {REC['fem2d']:.3e} / True)")
        print(f"[control] (Wc-FEM) shift P/2 miss {ids['wc']:.3e} over {ids['idx'].size} windows "
              f"> {D_BAND:.0e}  (W1:1972 {REC['wc']:.3e})")
        print(f"[printed] (Dc) m = 9 dropped {ids['dc']:.6f}  (W1:1973 {REC['dc']:.6f})")
        for m, h in enumerate(hs):
            print(f"[record] H_{m} <theta>(t*) = {h[-1]:.12e} K  (W1:{1957 + m} "
                  f"{REC['h_tstar'][m]:.12e})")
        print(f"[record] duty <theta> n=32 {duty[32]:.9f} n=36 {duty[36]:.9f} n=40 "
              f"{duty[40]:.9f} K  (W1:1967 "
              + " / ".join(f"{x:.9f}" for x in REC["duty_32_36_40"]) + ")")
        print(f"[record] duty max six-minute mean {ids['a_mx']:.9e} K at t_max = "
              f"{ids['a_tmx']:.4f} s; H_0 {ids['h0_mx']:.9e} K at {ids['h0_tmx']:.4f} s  "
              f"(W1:1974 {REC['a_mx']:.9e} / {REC['h0_mx']:.9e} at {REC['t_max']:.4f})")
        print(f"[anchor] (Cm) volume mean of the per-vertex six-minute field vs the windowed "
              f"series {cm:.3e} <= {D_BAND:.0e}  (predicted <= 1e-14)")
        print(f"[control] (Cm-c) window W - DT = {W - DT:.4f} s: miss {cmc:.3e} > {D_BAND:.0e}  "
              f"(predicted [5e-3, 6e-2])")
        print(f"[printed] source six-minute on-fraction max {np.nanmax(frac):.6f} at n = {n_frac} "
              f"(window [{t[n_frac] - W:.1f}, {t[n_frac]:.1f}] s) = {np.nanmax(frac) / 0.5:.4f} x D; "
              f"at n = 40 {frac[40]:.6f}  (predicted 0.651667 / 1.3033 / 0.348333)")
        print(f"[printed] volume-mean ratio a_max / H_0 max = {ratio:.4f}  (predicted "
              f"[0.52, 0.53], ~D); chi ratio {np.nanmax(frac) / 0.5:.4f}")
        print(f"[printed] windows: {int(ok.sum())} (n = {int(np.nonzero(ok)[0][0])} ... "
              f"{int(np.nonzero(ok)[0][-1])}), W / DT = {W / DT:.4f}")
        print(f"[printed] per-rank (owned P1 dofs, owned+ghost columns): {owned}")
        print(f"[printed] per-vertex: max six-minute mean {mx_max:.9e} K vs max theta(t*) "
              f"{th_tstar_max:.9e} K over {n_tot} owned vertices")
        print(f"[printed] vertices whose six-minute mean peaks before t*: {n_early} of {n_tot} "
              f"({n_early / n_tot:.4f}); at n = 36 (t = {t36:.1f} s): {n_36}; other early: "
              f"{n_early - n_36}; at t*: {n_tstar}  (predicted fraction [0, 0.5]; a zero is a "
              f"reading)")
        if partb is None:
            print("[Part B] input absent -- skipped (logs/mat15_s5/ not regenerated here)")
        else:
            dose, quad, b_mx, b_tmx, dmin = partb
            print(f"[Part B] 10 g ball-mean CEM43 per accepted W (ruling 361-1; not a safety "
                  f"statement): cem43 {dose:.15e} s, quad {quad:.15e} s  (W2:52 "
                  f"{REC['b_cem43_s']:.15e} / {REC['b_quad_s']:.15e})")
            print(f"[Part B] six-minute peak theta {b_mx:.15e} K at t_max = {b_tmx:.6f} s  "
                  f"(W2:57 {REC['b_peak_k']:.15e} at 1254.000000)")
            print(f"[Part B] dose {dmin:.6e} min per accepted W (ruling 361-1; not a safety "
                  f"statement)  (W2:62 {REC['b_dose_min']:.6e})", flush=True)

    # ---- asserts (every rank) ------------------------------------------------------
    assert ids["h0"] <= L_BAND, f"(H0) {ids['h0']:.3e}"
    assert ids["fem"] <= D_BAND, f"(FEM) {ids['fem']:.3e}"
    assert ids["fem2d"] <= D_BAND and ids["fem2d_eq"] is True, f"(FEM-2d) {ids['fem2d']:.3e}"
    assert ids["wc"] > D_BAND, f"(Wc-FEM) did not miss: {ids['wc']:.3e}"
    assert cm <= D_BAND, f"(Cm) {cm:.3e}"
    assert cmc > D_BAND, f"(Cm-c) did not miss: {cmc:.3e}"

    # ---- outputs ----------------------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)

    def _field(name, arr):
        f = fem.Function(V1, name=name)
        f.x.array[:] = arr
        f.x.scatter_forward()
        return f

    n_ok = [int(n) for n in np.nonzero(ok)[0]]
    steps = [(n * DT, {"theta_K": _field("theta_K", F[n]),
                       "theta_six_min_mean_K": _field("theta_six_min_mean_K", means_F[n])})
             for n in n_ok]
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_series_combined", sub, None, steps, comm=comm
    )
    q0.name = "q_h_k0"  # the XDMF dataset name is the Function's, not the dict key
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None,
        {"theta_duty_tstar_K": _field("theta_duty_tstar_K", F[-1]),
         "theta_six_min_mean_tstar_K": _field("theta_six_min_mean_tstar_K", means_F[-1]),
         "theta_six_min_max_K": _field("theta_six_min_max_K", mx_F),
         "theta_six_min_t_max_s": _field("theta_six_min_t_max_s", tmx_F),
         "q_h_k0": q0},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_window_history.csv"

        def _f(x):
            return "" if (x is None or (isinstance(x, float) and math.isnan(x))) else repr(float(x))

        res_full = np.full(t.size, np.nan)
        res_full[ok] = ids["res"]
        with open(csv_path, "w") as fh:
            fh.write("n,t_s,chi_step,chi_six_min_on_fraction,mean_theta_duty_K,mean_theta_H0_K,"
                     "six_min_duty_K,six_min_sum_Hm_K,rel_residual,six_min_H0_K\n")
            for n in range(t.size):
                chi = "" if n == 0 else repr(float(m17s3._chi_duty(t[n] - DT / 2.0)))
                fh.write(",".join([str(n), repr(float(t[n])), chi, _f(frac[n]), _f(duty[n]),
                                   _f(hs[0][n]), _f(a[n]), _f(ids["b"][n]), _f(res_full[n]),
                                   _f(ids["ra"][0][0][n])]) + "\n")
        adopt_host_ownership(OUTPUT_DIR, comm=MPI.COMM_SELF)
        print(f"\n[output] {series_path.name}: P1 theta_K and theta_six_min_mean_K at "
              f"{len(steps)} frames, t_n = n*{DT:.4f} s, n = {n_ok[0]} ... {n_ok[-1]}")
        print(f"[output] {xdmf_path.name}: P1 theta_duty_tstar_K, theta_six_min_mean_tstar_K, "
              f"theta_six_min_max_K, theta_six_min_t_max_s and DG0 q_h_k0 (t* = {T_STAR:.1f} s)")
        print(f"[output] {csv_path.name}: {t.size} rows")
        print("EX-110 anchors and controls: PASS")


if __name__ == "__main__":
    main()
