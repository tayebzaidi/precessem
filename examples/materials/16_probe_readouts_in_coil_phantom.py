"""Example (`EX-100`): what a probe log of the coil-heated gel phantom reports.

`mat:9` heats the coil-driven F-small gel phantom (`MAT-15` P1: `ANS-2`'s untuned
10 MHz loaded birdcage, single-leg drive k = 0, G = 1 W deposited, P1 thermal on
the 58 866-cell tag-3 submesh, CN + Rannacher, N = 40 to t* = 1254 s) but reads
theta only at t* and in 10 g balls. This one reads it the way a fibre-optic log
would: ball-averaged probe series (``MAT-19``'s virtual probe, an owned-cell
integral with one allreduce, never a point sample) at `ANS-9`'s SPEC probe
centres A-E, plus C' (the sigma_y image of C), the whole-phantom ball Omega and
the displaced control B+.

**Import, never restate.** The fixture is composed, in the order of
``MAT-19`` step 2's pytest fixture, from module-level names of
``tests/validation/test_mat19_probe_coil_phantom.py`` and the `MAT-15` module it
re-exports; every band and record is imported. Nothing in the gate modules was
edited.

Asserted (the step-2 identities, at imported bands): (S0) the imported objects
are the gate's, (A0) fixture identity and (M) containment, snapshot 0 reads 0.0,
(R0) the re-run is the gated run, (A) the whole-phantom probe is the volume mean,
(U4) bitwise rank agreement, (Sigma) the sigma_y mirror pair, (X1) the sweep
shares the series; controls (C1'') the displaced probe resolves at the width
band and (f'') the off-mirror pair; at -n 4 (U3c) the ghost-included Omega
misses (iff G_w > 0) and (W) the series agree with the -n 2 run at
``WIDTH_BAND`` (against that run's hash-guarded npz, never recomputed).
Printed only: the probe log, ball/point against the SPEC's row-4 points, the
placement sensitivity and an a-sweep at B.

Run it (two widths; the phantom is split across ranks only at -n 4)::

    ./scripts/run_examples.sh -e mat:16 -n 2 -t 590   # flagged run of record
    # then the same script at -n 4 (control window: asserts (W))

Outputs (``examples/materials/paraview_output/``):
``materials_16_probe_readouts_in_coil_phantom_combined.xdmf`` (41 time steps),
``materials_16_probe_readouts_in_coil_phantom_series.csv``,
``materials_16_probe_readouts_in_coil_phantom_probes.csv`` and
``materials_16_probe_readouts_in_coil_phantom_series.png``.
"""

from __future__ import annotations

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem
from dolfinx.mesh import meshtags

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_time_series,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import ProbeBall, probe_series  # noqa: E402
from fem_em_solver.thermal import probe as _probe_mod  # noqa: E402
from fem_em_solver.thermal.probe import _ball_average_on_cells  # noqa: E402

import tests.validation.test_mat15_coil_driven_heating_phantom as m15  # noqa: E402
import tests.validation.test_mat19_probe_coil_phantom as m19s2  # noqa: E402
import tests.validation.test_mat19_virtual_temperature_probe as m19  # noqa: E402
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    A,
    B_BAND,
    C4_COVARIANCE_BAND,
    E_BAND,
    F_CONTROL_FLOOR,
    N_MAIN,
    T_INF,
    T_STAR,
    V_BAND,
    V_SUB_RECORD_M3,
    _allsum,
    _run_heated,
    _transient,
)
from tests.validation.test_mat19_probe_coil_phantom import (  # noqa: E402
    A_PROBE,
    DT,
    N_SUB,
    NPZ4,
    ROW2_MEAN_K,
    ROW4_POINT_K,
    _probe_set,
    _sha,
)
from tests.validation.test_mat19_virtual_temperature_probe import WIDTH_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_16_probe_readouts_in_coil_phantom"
NPZ_PATH = "/workspace/logs/ex100/mat16-n2.npz"
A_SWEEP = (0.5 * A_PROBE, A_PROBE)  # example parameters at B; 5 mm is the gate's A_PROBE
SPEC_TIMES = (10, 20, 30, 40)
SAFE = {"C'": "Cprime", "B+": "Bplus"}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _cell_centroids(sub):
    return sub.geometry.x[sub.geometry.dofmap].mean(axis=1)


def _inside(sub, centre, a):
    return np.sum((_cell_centroids(sub) - np.asarray(centre, dtype=float)) ** 2, axis=1) <= a * a


def _ball_mask_dg0(sub, centre, a, name):
    """DG0 1/0 mask of cells whose centroid lies inside the ball (float64 numpy)."""
    V0 = fem.functionspace(sub, ("DG", 0))
    f = fem.Function(V0, name=name)
    ncell = sub.geometry.dofmap.shape[0]
    dofs = np.array([V0.dofmap.cell_dofs(c)[0] for c in range(ncell)], dtype=np.int64)
    f.x.array[dofs] = np.where(_inside(sub, centre, a), 1.0, 0.0)
    f.x.scatter_forward()
    return f


def _figure_tags(sub, centres):
    cen = _cell_centroids(sub)
    vals = np.ones(cen.shape[0], dtype=np.int32)
    for tag, k in zip((11, 12, 13, 14, 15), "ABCDE"):
        ins = _inside(sub, centres[k], A_PROBE)
        vals[ins & (vals == 1)] = tag
    idx = np.arange(cen.shape[0], dtype=np.int32)
    return meshtags(sub, sub.topology.dim, idx, vals)


def main() -> None:  # noqa: C901
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )
    evaluate_w = comm.size in (2, 4)

    _say(comm, "=" * 78)
    _say(comm, "EX-100 -- probe readouts at the benchmark's centres, coil-heated gel phantom")
    _say(comm, "=" * 78)

    # ---- (S0) structural --------------------------------------------------------
    s0 = {
        "B_BAND": m19s2.B_BAND is m15.B_BAND,
        "WIDTH_BAND": m19s2.WIDTH_BAND is m19.WIDTH_BAND,
        "_run_heated": m19s2._run_heated is m15._run_heated,
        "C4_COVARIANCE_BAND": m19s2.C4_COVARIANCE_BAND is m15.C4_COVARIANCE_BAND,
    }
    _say(comm, f"[S0] imported objects identical to the gates': {s0}")
    assert all(s0.values()), s0
    if not evaluate_w:
        _say(comm, f"[note] -n {comm.size}: readings printed, (W) not evaluated (run at -n 2 then -n 4)")

    # ---- fixture, in step 2's order -------------------------------------------
    t0 = time.perf_counter()
    out = _run_heated(return_fields=True)
    t_fix = time.perf_counter() - t0
    sub, ft = out["sub"], out["ft"]
    q_h = out["q_h"][0]
    v_sub = out["v_sub"]
    n_sub = sub.topology.index_map(3).size_global
    names, order, probes, centres = _probe_set()
    j = {k: order.index(k) for k in order}

    write_setup_figure(
        sub,
        _figure_tags(sub, centres),
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: "gel phantom (tag 3 submesh, G = 1 W)",
            **{t: f"probe {k} (a = {A_PROBE * 1e3:.0f} mm) at ({centres[k][0] * 1e3:+.0f}, "
                  f"{centres[k][1] * 1e3:+.0f}, {centres[k][2] * 1e3:+.0f}) mm"
               for t, k in zip((11, 12, 13, 14, 15), "ABCDE")},
        },
        translucent_tags=(1,),
        hide_tags=(),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        window_size=(1400, 630),
        title="mat:16 -- ball-averaged probes at ANS-9's SPEC centres in the coil-heated gel phantom (MAT-19 step 2)",
        comm=comm,
    )

    t1 = time.perf_counter()
    res = _transient(sub, ft, q_h, dt=DT, n_steps=N_MAIN, theta=0.5, startup_half_steps=2,
                     snapshot_steps=range(N_MAIN + 1))
    d = res.snapshots[N_MAIN] - out["T_star"][0]
    th_star = out["T_star"][0] - T_INF
    r0 = float(np.sqrt(np.real(_allsum(comm, ufl.inner(d, d) * ufl.dx(domain=sub))))
               / np.sqrt(np.real(_allsum(comm, ufl.inner(th_star, th_star) * ufl.dx(domain=sub)))))
    V = res.snapshots[N_MAIN].function_space
    th = {}
    for n in range(N_MAIN + 1):
        f = fem.Function(V)
        f.x.array[:] = res.snapshots[n].x.array - T_INF
        f.x.scatter_forward()
        th[n] = f
    t_rerun = time.perf_counter() - t1

    t2 = time.perf_counter()
    s = probe_series(th, DT, probes)
    v = np.real(s.values)
    sweep = [ProbeBall(centres["B"], a) for a in A_SWEEP]
    s_x = probe_series(th, DT, sweep)
    vx = np.real(s_x.values)
    t_series = time.perf_counter() - t2
    means = np.array([float(np.real(_allsum(comm, th[n] * ufl.dx(domain=sub)))) / v_sub
                      for n in range(N_MAIN + 1)])

    imap = sub.topology.index_map(3)
    part = comm.allgather((int(imap.size_local), int(imap.num_ghosts)))
    o_w = sum(1 for so, _g in part if so > 0)
    g_w = sum(g for _so, g in part)
    ghost = None
    if comm.size > 1:
        allc = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
        gval, gvol = _ball_average_on_cells(th[N_MAIN], probes[j["Omega"]], allc)
        ghost = (abs(gvol / v_sub - 1.0), abs(gval / means[N_MAIN] - 1.0))
    allrows = comm.allgather(s.values[N_MAIN].copy())
    u4 = all(np.array_equal(allrows[0], r) for r in allrows)

    # ---- (A0) -----------------------------------------------------------------------
    err_V = abs(v_sub / V_SUB_RECORD_M3 - 1.0)
    m_ok = all(
        np.hypot(c[0], c[1]) + a <= A and abs(c[2]) + a <= PHANTOM_HEIGHT / 2.0
        for k in order if k != "Omega"
        for c, a in [(centres[k], A_PROBE)] + [(centres["B"], x) for x in A_SWEEP] * (k == "B")
    )
    _say(comm, f"\n[A0] drive by_k[0] = {out['by_k'][0]}; P_1 = {out['P'][0]:.12e} W; "
               f"submesh cells {n_sub} (record {N_SUB}); V_sub rel {err_V:.3e} (bound {V_BAND:.0e}); "
               f"(M) every ball inside the phantom: {m_ok}; width -n {comm.size}")
    _say(comm, f"[part] -n {comm.size}: sub owned {[p[0] for p in part]} ghosts {[p[1] for p in part]} "
               f"O_w {o_w} G_w {g_w} (D4 record at -n 4: owned 23401 / 35465 / 0 / 0, ghosts 701 / 677 / 0 / 0)")
    assert out["by_k"][0] == "P1"
    assert n_sub == N_SUB
    assert err_V <= V_BAND
    assert m_ok
    assert np.all(s.values[0] == 0.0), "snapshot 0 must read exactly 0.0"
    _say(comm, f"[R0] |theta40_rerun - theta*_gated| / |theta*_gated| = {r0:.3e} (bound {E_BAND:.0e}; "
               f"record 2.531e-15)")
    assert r0 <= E_BAND

    # ---- (A) ------------------------------------------------------------------------
    PO = v[:, j["Omega"]]
    relA = float(np.max(np.abs(PO[1:] - means[1:]) / np.abs(means[1:])))
    volA = abs(s.volumes[j["Omega"]] / v_sub - 1.0)
    _say(comm, f"[A] P_Omega vs int theta/V_sub: max over n = 1..40 {relA:.3e}; V_Omega/V_sub - 1 = "
               f"{volA:.3e} (bound {B_BAND:.0e}; record 2.190e-14 at -n 2, 6.825e-15 at -n 4)")
    _say(comm, "[A] P_Omega at the SPEC times vs row 2: " + "  ".join(
        f"n={n} {PO[n]:.6f} (row2 {ROW2_MEAN_K[n]:.6f})" for n in SPEC_TIMES))
    assert relA <= B_BAND and volA <= B_BAND
    _say(comm, f"[U4] bitwise rank agreement at n = 40: {u4} (-n {comm.size})")
    if comm.size > 1:
        assert u4

    # ---- (Sigma), (C1''), (f'') -----------------------------------------------------
    sig = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["C'"]]) / v[N_MAIN, j["C"]]
    c1 = abs(v[N_MAIN, j["B+"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    ff = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(comm, f"[Sigma] |C - C'|/C at t* = {sig:.4e} (bound {C4_COVARIANCE_BAND:.0e}; record 5.9190e-04)  ASSERTED")
    _say(comm, f"[control] (C1'') |B+ - B|/B at t* = {c1:.4e} (asserted > {WIDTH_BAND:.0e}; record 4.3978e-03)")
    _say(comm, f"[control] (f'') |C - B|/B at t* = {ff:.4f} (asserted > {F_CONTROL_FLOOR:.0e}; record 0.8018)")
    assert sig <= C4_COVARIANCE_BAND
    assert c1 > WIDTH_BAND
    assert ff > F_CONTROL_FLOOR

    # ---- (X1) -----------------------------------------------------------------------
    x1 = float(np.max(np.abs(vx[1:, 1] / v[1:, j["B"]] - 1.0)))
    _say(comm, f"[X1] sweep a = {A_SWEEP[1] * 1e3:.1f} mm column vs main B series: max rel {x1:.3e} "
               f"(bound {B_BAND:.0e}); steps == range({N_MAIN + 1}): {list(s_x.steps) == list(range(N_MAIN + 1))}")
    assert x1 <= B_BAND
    assert list(s_x.steps) == list(range(N_MAIN + 1))

    # ---- (U3c) ----------------------------------------------------------------------
    if ghost is not None:
        if g_w > 0:
            _say(comm, f"[control] (U3c) ghost-included Omega: volume miss {ghost[0]:.4e} (asserted > "
                       f"{B_BAND:.0e}; record 2.6047e-02 at -n 4), value miss {ghost[1]:.4e}")
            assert ghost[0] > B_BAND
        else:
            _say(comm, f"[control] (U3c) inapplicable at -n {comm.size}: G_w = 0 (volume miss {ghost[0]:.4e} printed)")

    # ---- printed: probe log, (P4), (P5), (P-a) ----------------------------------------
    _say(comm, "\n[log] probe log, theta_bar [K] per 1 W absorbed (T = T_inf + theta_bar; linear in power); "
               "t = " + " / ".join(f"{n * DT:.1f} s" for n in SPEC_TIMES))
    for k in names + ["C'"]:
        _say(comm, f"[log]   {k:2s} theta " + " ".join(f"{v[n, j[k]]:.7f}" for n in SPEC_TIMES)
                   + "   T " + " ".join(f"{T_INF + v[n, j[k]]:.5f}" for n in SPEC_TIMES))
    _say(comm, "[P4] ball(t*) / row-4 point - 1: " + "  ".join(
        f"{k} {v[N_MAIN, j[k]] / ROW4_POINT_K[k] - 1.0:+.4e}" for k in names)
        + "  (records A +7.1794e-03 B -1.2383e-02 C +1.8067e-02 D -4.3455e-03 E +2.8398e-03; the SPEC's "
          "row 4 accepts a ball of radius <= 5 mm or the point; printed only)")
    _say(comm, f"[P5] placement: 0.1 mm radial shift of B moves the reading by {c1:.4e} "
               f"(about {c1 / 0.1 * 100:.2f} %/mm); a sensitivity, not a gate")
    _say(comm, "[P-a] probe volume at B: a [mm], V_h/(4 pi a^3/3) - 1, theta_bar(t*)/point - 1 "
               "(not gated; a > 5 mm would exceed the SPEC allowance)")
    vv = s_x.volumes / s_x.nominal_volumes - 1.0
    for i, a in enumerate(A_SWEEP):
        _say(comm, f"[P-a]   {a * 1e3:4.1f}  {vv[i]:+.3e}  {vx[N_MAIN, i] / ROW4_POINT_K['B'] - 1.0:+.4e}")
    _say(comm, "[vol] ball V_h/(4 pi a^3/3) - 1: " + "  ".join(
        f"{k} {s.volumes[j[k]] / s.nominal_volumes[j[k]] - 1.0:+.3e}" for k in order if k != "Omega"))
    _say(comm, f"[imag] max_imag_ratio {s.max_imag_ratio:.3e} (predicted 0)")

    # ---- the npz and (W) ---------------------------------------------------------------
    me = Path(__file__).resolve()
    hashes = (_sha(me), _sha(m19s2.__file__), _sha(_probe_mod.__file__))
    if comm.size == 2:
        if comm.rank == 0:
            os.makedirs(os.path.dirname(NPZ_PATH), exist_ok=True)
            np.savez(NPZ_PATH, series=v, volumes=s.volumes, cells=np.array([n_sub]),
                     sha_script=hashes[0], sha_module=hashes[1], sha_probe=hashes[2])
            print(f"[W] wrote the -n 2 npz (sha script {hashes[0][:12]}, module {hashes[1][:12]}, "
                  f"probe {hashes[2][:12]})", flush=True)
        comm.Barrier()
        if os.path.exists(NPZ4):
            z4 = np.load(NPZ4)
            if str(z4["sha_module"]) == hashes[1] and str(z4["sha_probe"]) == hashes[2]:
                _say(comm, "[P-W'] this run vs the gate's -n 4 series: max rel "
                           f"{float(np.max(np.abs(v[1:] / z4['series'][1:] - 1.0))):.3e} (predicted <= 1e-10; printed)")
            else:
                _say(comm, "[P-W'] gate npz absent or stale -- not read")
        else:
            _say(comm, "[P-W'] gate npz absent or stale -- not read")
    elif comm.size == 4:
        if not os.path.exists(NPZ_PATH):
            raise RuntimeError("(W) the -n 2 npz is missing: run the flagged -n 2 run first; never recompute")
        z = np.load(NPZ_PATH)
        if (str(z["sha_script"]), str(z["sha_module"]), str(z["sha_probe"])) != hashes:
            raise RuntimeError("(W) npz hashes differ from this script / the step-2 module / thermal/probe.py: "
                               "re-run the -n 2 run")
        cells_eq = int(z["cells"][0]) == n_sub
        rel = np.abs(v[1:] / z["series"][1:] - 1.0)
        w = float(np.max(rel))
        _say(comm, f"[W] -n 4 vs -n 2: cells equal {cells_eq}; series max rel {w:.3e} (band {WIDTH_BAND:.0e}; "
                   f"predicted <= 1e-10, step-2 record 3.714e-11); per probe max: " + "  ".join(
                       f"{k} {np.max(rel[:, j[k]]):.2e}" for k in order))
        assert cells_eq
        assert w <= WIDTH_BAND

    # ---- XDMF ------------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    masks = {}
    for k in names + ["C'", "B+"]:
        nm = f"probe_{SAFE.get(k, k)}_mask"
        masks[nm] = _ball_mask_dg0(sub, centres[k], A_PROBE, nm)
    q_h.name = "q_h"
    frames = []
    for n in range(N_MAIN + 1):
        th[n].name = "theta"
        frames.append((n * DT, {"theta": th[n], "q_h": q_h, **masks}))
    xdmf_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, frames, comm=comm,
    )

    series_csv = OUTPUT_DIR / f"{BASENAME}_series.csv"
    probes_csv = OUTPUT_DIR / f"{BASENAME}_probes.csv"
    png = OUTPUT_DIR / f"{BASENAME}_series.png"
    if comm.rank == 0:
        with series_csv.open("w", newline="") as fh:
            w_ = csv.writer(fh)
            head = ["n", "t_s"]
            for k in order:
                head += [f"theta_{SAFE.get(k, k)}_K", f"T_{SAFE.get(k, k)}_K"]
            w_.writerow(head)
            for n in range(N_MAIN + 1):
                row = [n, f"{n * DT:.6f}"]
                for k in order:
                    row += [f"{v[n, j[k]]:.9e}", f"{T_INF + v[n, j[k]]:.9e}"]
                w_.writerow(row)
        with probes_csv.open("w", newline="") as fh:
            w_ = csv.writer(fh)
            w_.writerow(["label", "cx_m", "cy_m", "cz_m", "a_m", "V_h_over_V_minus_1", "ball_over_point_minus_1_tstar"])
            for k in order:
                a_k = m19s2.OMEGA_RADIUS if k == "Omega" else A_PROBE
                r4 = (f"{v[N_MAIN, j[k]] / ROW4_POINT_K[k] - 1.0:.9e}" if k in ROW4_POINT_K else "")
                w_.writerow([SAFE.get(k, k), *(f"{x:.9e}" for x in centres[k]), f"{a_k:.3e}",
                             f"{s.volumes[j[k]] / s.nominal_volumes[j[k]] - 1.0:.9e}", r4])
            for i, a in enumerate(A_SWEEP):
                w_.writerow([f"B_a{a * 1e3:.1f}mm", *(f"{x:.9e}" for x in centres["B"]), f"{a:.3e}",
                             f"{vv[i]:.9e}", f"{vx[N_MAIN, i] / ROW4_POINT_K['B'] - 1.0:.9e}"])

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ts = np.arange(N_MAIN + 1) * DT
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for k in names:
            (line,) = ax.plot(ts, v[:, j[k]], "-", label=f"probe {k}")
            ax.plot([T_STAR], [ROW4_POINT_K[k]], "o", color=line.get_color(), ms=6)
        ax.set_xlabel("t [s]")
        ax.set_ylabel("theta_bar [K per W absorbed]")
        ax.set_title("ball-averaged probes (markers: SPEC row-4 point records at t*)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(png, dpi=110)
        plt.close(fig)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf_path} -- 41 steps: theta, q_h, probe masks")
    _say(comm, f"[output] {series_csv}")
    _say(comm, f"[output] {probes_csv}")
    _say(comm, f"[output] {png}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size} (fixture "
               f"{t_fix:.1f} s, rerun {t_rerun:.1f} s, series {t_series:.1f} s)")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
