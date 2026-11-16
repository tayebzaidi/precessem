"""Example (`EX-90`): what a temperature probe in a heated gel ball reports, and
why its stated volume matters.

`mat:4` stops at the whole-ball mean and the interior L2 norm; `mat:11` prints
a centre history. This one reads the heated gel ball through ``MAT-19``'s
virtual temperature probe: a ball-averaged temperature against time, computed
as an integral over owned cells with one allreduce, never as a point sample.
One EM solve on `MAT-4`'s lossy sphere (R = 10 mm, sigma 0.57 S/m, 64 MHz,
imposed field) heats the 19 496-cell sphere submesh (Robin h = 100 W/m^2/K,
CN + Rannacher, N = 40 to t* = 139.33 s).

* Three probes of stated radius a = 2 mm (centres 0, 3.5 mm on x, 6.5 mm
  off-axis) are set beside the exact ball mean of the Robin eigen-series
  (asserted by the gate to 1 %) and beside the point value a bare comparison
  would use.
* At the centre the same run is read again with a = 1, 4 and 8 mm: the volume
  offset (ball - point) grows as a^2 and is invisible to a 1 % band at 2 mm,
  but not at 8 mm.

**Import, never restate** (the `ANS-1` rule). The fixture is
``_build_run()`` and the gate body is ``_check_probe_series(run)``, both lifted
(rule (a), additive return keys only) from
``tests/validation/test_mat19_virtual_temperature_probe.py``; the helpers
``_ball_mean``, ``_point``, ``_q_cf`` and every constant and band come from
that module or from the `MAT-10` / `MAT-12` modules it re-exports.

Asserted: the gate's own asserts inside the lifted body ((M), (A0), snapshot 0,
(T) <= ``P_BAND``, (C1) > ``P_BAND``), re-asserted here on the returned values;
(X0) sweep containment |c0| + a <= 0.9 R; (X1) the sweep call's first column
equals the gate series. Printed, not asserted: the records, the sweep, the
volume effect, and (C3), the volume-omitted comparand at a = 8 mm.

Run it::

    ./run_examples.sh -e mat:13 -n 2 -t 300

Outputs (``examples/materials/paraview_output/``):
``materials_13_virtual_probe_in_gel_ball_combined.xdmf`` (41 time steps),
``materials_13_virtual_probe_in_gel_ball_series.csv``,
``materials_13_virtual_probe_in_gel_ball_probes.csv`` and
``materials_13_virtual_probe_in_gel_ball_series.png``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
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

import tests.validation.test_mat10_sar_to_heat_sphere as mat10  # noqa: E402
import tests.validation.test_mat12_transient_sar_heating_sphere as mat12  # noqa: E402
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND  # noqa: E402
from tests.validation.test_mat19_virtual_temperature_probe import (  # noqa: E402
    A_PROBE,
    CENTRES,
    DISPLACED,
    DT,
    H,
    K,
    N_MAIN,
    N_SUB,
    P_BAND,
    R,
    T_STAR,
    _ball_mean,
    _build_run,
    _check_probe_series,
    _point,
    _q_cf,
    _theta_ss0,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_13_virtual_probe_in_gel_ball"
A_SWEEP = (1.0e-3, 4.0e-3, 8.0e-3)  # example parameters at c0, not records
W1 = "20261005T095808Z_MAT-19-step1-W1.log"

#: Gate records, `20261005T095808Z_MAT-19-step1-W1.log` (line in the comment).
#: Printed beside this run, a difference reported and never asserted.
RECORD = {
    "n_sub": 19496,                                      # W1:161
    "vol": (6.289e-04, 4.055e-04, 5.064e-04, 1.569e-04),  # W1:161
    "P_t": 1.6845e-03,                                   # W1:160
    "worstT": 1.7315e-03,                                # W1:203
    "Pstar": ("5.918962e-02", "5.681504e-02", "5.020215e-02"),  # W1:202
    "qh": (3.7300e-02, 3.6500e-02, 3.4656e-02),          # W1:207
    "missC1": 5.4207e-02,                                # W1:209
    "ownrel": 1.1349e-03,                                # W1:209
    "c2": (-6.2905e-03, -8.1485e-03, -1.3531e-02),       # W1:210
    "c2n20": -1.1388e-02,                                # W1:210
}
#: Predictions written before the run (review 2026-10-05 10:30), not records.
PRED = {
    "off_exact_tstar": {1.0e-3: -0.180e-2, 2.0e-3: -0.726e-2, 4.0e-3: -2.994e-2, 8.0e-3: -13.368e-2},
    "off_exact_n20": {1.0e-3: -0.079e-2, 2.0e-3: -0.324e-2, 4.0e-3: -1.450e-2, 8.0e-3: -8.444e-2},
    "off_fem_range": {1.0e-3: (-0.30e-2, 0.05e-2), 4.0e-3: (-3.3e-2, -2.6e-2),
                      8.0e-3: (-13.8e-2, -12.9e-2)},
    "ratio_exact": {1.0e-3: 0.248, 4.0e-3: 4.12, 8.0e-3: 18.4},
    "c3": (12.9e-2, 13.8e-2),
    "own_bar": 3.0e-3,
    "vol_abs": {1.0e-3: 5.0e-3, 4.0e-3: 2.0e-3, 8.0e-3: 2.0e-3},
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _cell_centroids(sub):
    return sub.geometry.x[sub.geometry.dofmap].mean(axis=1)


def _ball_mask_dg0(sub, centre, a, name):
    """DG0 1/0 mask of cells whose centroid lies inside the ball (float64 numpy)."""
    V0 = fem.functionspace(sub, ("DG", 0))
    f = fem.Function(V0, name=name)
    cen = _cell_centroids(sub)
    inside = np.sum((cen - np.asarray(centre, dtype=float)) ** 2, axis=1) <= a * a
    ncell = cen.shape[0]
    dofs = np.array([V0.dofmap.cell_dofs(c)[0] for c in range(ncell)], dtype=np.int64)
    f.x.array[dofs] = np.where(inside, 1.0, 0.0)
    f.x.scatter_forward()
    return f


def _figure_tags(sub):
    cen = _cell_centroids(sub)
    vals = np.ones(cen.shape[0], dtype=np.int32)
    for tag, c in ((11, CENTRES[0]), (12, CENTRES[1]), (13, CENTRES[2])):
        inside = np.sum((cen - np.asarray(c, dtype=float)) ** 2, axis=1) <= A_PROBE**2
        vals[inside & (vals == 1)] = tag
    idx = np.arange(cen.shape[0], dtype=np.int32)
    return meshtags(sub, sub.topology.dim, idx, vals)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )
    if comm.size != 2:
        raise RuntimeError("run at -n 2: the lifted gate body reads/writes the -n 2 npz")

    _say(comm, "=" * 78)
    _say(comm, "EX-90 -- a virtual temperature probe in a heated gel ball")
    _say(comm, "=" * 78)

    same = {
        "P_BAND": P_BAND is mat10.P_BAND, "B_BAND": B_BAND is mat10.B_BAND,
        "R": R is mat10.R, "N_MAIN": N_MAIN is mat12.N_MAIN, "T_STAR": T_STAR is mat12.T_STAR,
    }
    _say(comm, f"[gate] imported bands/constants identical to MAT-10/MAT-12 objects: {same}; "
               f"P_BAND {P_BAND:.0e}, B_BAND {B_BAND:.0e}, R {R}, N {N_MAIN}, t* {T_STAR:.4f} s")
    assert all(same.values()), same

    # ---- the gate fixture, once ------------------------------------------------
    run = _build_run()
    sub = run["sub"]
    q = run["q_cf"]
    assert abs(q - _q_cf()) == 0.0

    write_setup_figure(
        sub,
        _figure_tags(sub),
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: "gel sphere (sigma 0.57 S/m, Robin h = 100)",
            11: "probe c0 (a = 2 mm)",
            12: "probe c1 (a = 2 mm, 3.5 mm on x)",
            13: "probe c2 (a = 2 mm, 6.5 mm off-axis)",
        },
        translucent_tags=(1,),
        hide_tags=(),
        slice_normal=tuple(np.array([0.0, -1.0, 1.0]) / np.sqrt(2.0)),
        slice_origin=(0.0, 0.0, 0.0),
        window_size=(1400, 630),
        title="mat:13 -- ball-averaged probes in the heated gel ball (MAT-19 on MAT-12's fixture)",
        comm=comm,
    )

    # ---- the gate body, through the lift --------------------------------------
    chk = _check_probe_series(run)
    s_cf = run["s_cf"]

    # ---- anchors on the returned values, beside the W1 records --------------------
    _say(comm, "\n[anchor] re-asserted on the lifted gate's returned values; records "
               f"from {W1}")
    vol = s_cf.volumes / s_cf.nominal_volumes - 1.0
    _say(comm, f"[anchor] (A0) submesh cells {run['n_sub']} (record {RECORD['n_sub']}); V_h/V-1 "
               + " ".join(f"{v:+.3e}" for v in vol) + " (record "
               + " ".join(f"{v:+.3e}" for v in RECORD["vol"]) + ")")
    assert run["n_sub"] == N_SUB == RECORD["n_sub"]
    _say(comm, f"[anchor] (P_t) {run['err_P']:.4e} (record {RECORD['P_t']:.4e}, diff "
               f"{run['err_P'] - RECORD['P_t']:+.1e})")
    P = chk["P"]
    pstar = tuple(f"{P[N_MAIN, j]:.6e}" for j in range(3))
    _say(comm, f"[anchor] (T) readings at t*: {' '.join(pstar)} (record {' '.join(RECORD['Pstar'])}; "
               f"{'identical' if pstar == RECORD['Pstar'] else 'DIFFER'})")
    _say(comm, f"[anchor] (T) max over 3 probes x 40 times {chk['worstT']:.4e} (record "
               f"{RECORD['worstT']:.4e}, diff {chk['worstT'] - RECORD['worstT']:+.1e}); "
               f"ASSERTED <= {P_BAND:.0e}")
    assert chk["worstT"] <= P_BAND
    c1 = float(np.max(chk["missC1"]))
    _say(comm, f"[anchor] (C1) displaced probe vs nominal theta_bar(6.5 mm): max miss {c1:.4e} (record "
               f"{RECORD['missC1']:.4e}, diff {c1 - RECORD['missC1']:+.1e}); ASSERTED > {P_BAND:.0e}; "
               f"vs own theta_bar max {float(np.max(chk['ownrel'])):.4e} (record {RECORD['ownrel']:.4e})")
    assert c1 > P_BAND
    qh = [float(np.real(run["s_h"].values[N_MAIN, j]) / chk["ref"][N_MAIN, j] - 1.0) for j in range(3)]
    _say(comm, "[anchor] q_h series at t*: " + " ".join(f"{v:+.4e}" for v in qh) + " (record "
               + " ".join(f"{v:+.4e}" for v in RECORD["qh"]) + "); this is MAT-4's SAR offset, not the probe")
    _say(comm, "[anchor] (C2) ball series vs point value at t*: "
               + " ".join(f"{v:+.4e}" for v in chk["c2"]) + " (record "
               + " ".join(f"{v:+.4e}" for v in RECORD["c2"]) + f"); c2 at n=20 {chk['c2n20']:+.4e} "
               f"(record {RECORD['c2n20']:+.4e})")

    # ---- (X0), (X1): the sweep at the centre ------------------------------------
    d0 = float(np.linalg.norm(CENTRES[0]))
    for a in A_SWEEP:
        assert d0 + a <= 0.9 * R, f"(X0) a = {a} leaves r < 0.9 R"
    _say(comm, f"\n[sweep] (X0) |c0| + a <= 0.9 R holds for a = {[a * 1e3 for a in A_SWEEP]} mm "
               f"(largest {d0 + max(A_SWEEP):.4f} <= {0.9 * R:.4f} m)  ASSERTED")
    sweep = [ProbeBall(CENTRES[0], A_PROBE)] + [ProbeBall(CENTRES[0], a) for a in A_SWEEP]
    s_x = probe_series(run["snap_cf"], DT, sweep)
    Px = np.real(s_x.values)
    x1 = float(np.max(np.abs(Px[1:, 0] / np.real(s_cf.values[1:, 0]) - 1.0)))
    _say(comm, f"[sweep] (X1) sweep column 0 vs gate series: max rel {x1:.3e} (predicted 0.0, "
               f"ASSERTED <= {B_BAND:.0e}); steps == range({N_MAIN + 1}): {list(s_x.steps) == list(range(N_MAIN + 1))}")
    assert x1 <= B_BAND
    assert list(s_x.steps) == list(range(N_MAIN + 1))

    ts = np.arange(N_MAIN + 1) * DT
    pt0 = np.array([float(_point(0.0, t, q)[0]) for t in ts])
    ex = {a: np.array([_ball_mean(0.0, a, t, q) for t in ts]) for a in (A_PROBE,) + A_SWEEP}
    sv = s_x.volumes / s_x.nominal_volumes - 1.0
    labels = (A_PROBE,) + A_SWEEP
    ss0 = _theta_ss0(q)
    _say(comm, "\n[sweep] a [mm]  V_h/V-1    max|P/theta_bar-1|  at t*      P/theta-1 at t* "
               "(exact; pred. FEM range)   n=20 exact   steady limit")
    off_fem, off_ex = {}, {}
    for j, a in enumerate(labels):
        own_max = float(np.max(np.abs(Px[1:, j] / ex[a][1:] - 1.0)))
        own_t = float(Px[N_MAIN, j] / ex[a][N_MAIN] - 1.0)
        off_fem[a] = float(Px[N_MAIN, j] / pt0[N_MAIN] - 1.0)
        off_ex[a] = float(ex[a][N_MAIN] / pt0[N_MAIN] - 1.0)
        n20 = float(ex[a][20] / pt0[20] - 1.0)
        steady = -q * a**2 / (10.0 * K) / ss0
        rng = PRED["off_fem_range"].get(a)
        rng_s = f"[{rng[0] * 100:+.2f}, {rng[1] * 100:+.2f}] %" if rng else "record -0.629 %"
        _say(comm, f"[sweep] {a * 1e3:5.1f}  {sv[j]:+.3e}  {own_max:.3e}  {own_t:+.3e}   "
                   f"{off_fem[a] * 100:+.3f} % (exact {off_ex[a] * 100:+.3f} % pred "
                   f"{PRED['off_exact_tstar'][a] * 100:+.3f} %; {rng_s})   "
                   f"{n20 * 100:+.3f} % (pred {PRED['off_exact_n20'][a] * 100:+.3f} %)   "
                   f"{steady * 100:+.3f} %")
        lim = PRED["own_bar"] if a != A_PROBE else P_BAND
        _say(comm, f"[sweep]        vs own theta_bar max {own_max:.3e} "
                   f"{'within' if own_max <= lim else 'OUTSIDE'} predicted {lim:.0e}; |V_h/V-1| "
                   f"{'within' if a == A_PROBE or abs(sv[j]) <= PRED['vol_abs'][a] else 'OUTSIDE'} prediction"
                   + ("" if a == A_PROBE else f" {PRED['vol_abs'][a]:.0e}"))
    _say(comm, "[sweep] NOT GATED: MAT-19's (T) covers a = 2 mm only; the sweep readings are printed. "
               "Below a ~ h_in (= R/10 = 1 mm) the volume offset is not resolved on this mesh.")
    for a in A_SWEEP:
        ratio_fem = off_fem[a] / off_fem[A_PROBE]
        ratio_ex = off_ex[a] / off_ex[A_PROBE]
        _say(comm, f"[a^2] offset ratio to a = 2 mm at t*, a = {a * 1e3:.0f} mm: FEM {ratio_fem:.3f}, "
                   f"exact {ratio_ex:.3f} (pred {PRED['ratio_exact'][a]}); steady part exactly "
                   f"{(a / A_PROBE) ** 2:.2f}")
    d = [float(np.linalg.norm(c)) for c in CENTRES]
    for j in (1, 2):
        e = _ball_mean(d[j], A_PROBE, T_STAR, q) / float(_point(d[j], T_STAR, q)[0]) - 1.0
        _say(comm, f"[c{j}] a = 2 mm exact theta_bar/theta - 1 at t*: {e:+.4e} (pred "
                   f"{'-0.865' if j == 1 else '-1.259'} %), FEM (C2) {chk['c2'][j]:+.4e}")

    # (C3): the volume-omitted comparand at a = 8 mm, predicted, printed
    P8 = Px[:, 3]
    c3_all = np.abs(P8[1:] / pt0[1:] - 1.0)
    c3 = float(np.max(c3_all))
    lo, hi = PRED["c3"]
    _say(comm, f"\n[control] (C3) a = 8 mm, ball probe vs point value theta(0, t_n): max over n "
               f"{c3 * 100:.3f} % at n = {int(np.argmax(c3_all)) + 1} (n=40: {c3_all[-1] * 100:.3f} %); "
               f"predicted {lo * 100:.1f}..{hi * 100:.1f} %, ceiling -17.45 %; band {P_BAND:.0e} -> "
               f"{c3 / P_BAND:.1f}x the band; {'MISSES' if c3 > P_BAND else 'DOES NOT MISS'} "
               f"(printed, not asserted)")
    if c3 <= P_BAND:
        for n in range(N_MAIN + 1):
            _say(comm, f"[control]   n={n:2d} P8 = {P8[n]:.6e}  theta = {pt0[n]:.6e}")

    # ---- XDMF: P1 theta_h / theta_exact / theta_error, DG0 masks ---------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = run["snap_cf"][0].function_space
    masks = {
        "probe_c0_mask": _ball_mask_dg0(sub, CENTRES[0], A_PROBE, "probe_c0_mask"),
        "probe_c1_mask": _ball_mask_dg0(sub, CENTRES[1], A_PROBE, "probe_c1_mask"),
        "probe_c2_mask": _ball_mask_dg0(sub, CENTRES[2], A_PROBE, "probe_c2_mask"),
        "probe_c2_displaced_mask": _ball_mask_dg0(sub, DISPLACED, A_PROBE, "probe_c2_displaced_mask"),
        "probe_c0_a8mm_mask": _ball_mask_dg0(sub, CENTRES[0], A_SWEEP[-1], "probe_c0_a8mm_mask"),
    }
    frames = []
    for n in range(N_MAIN + 1):
        th = run["snap_cf"][n]
        th.name = "theta_h"
        ex_f = fem.Function(V1, name="theta_exact")

        def ev(x, t=n * DT):
            r = np.sqrt(np.sum(x**2, axis=0))
            out = np.empty(r.shape)
            for s in range(0, r.size, 20000):
                out[s:s + 20000] = _point(r[s:s + 20000], t, q)
            return out

        ex_f.interpolate(ev)
        ex_f.x.scatter_forward()
        er = fem.Function(V1, name="theta_error")
        er.x.array[:] = th.x.array - ex_f.x.array
        er.x.scatter_forward()
        frames.append((n * DT, {"theta_h": th, "theta_exact": ex_f, "theta_error": er, **masks}))
    xdmf_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, frames, comm=comm,
    )

    # ---- tables and figure (rank 0, no collective inside) ---------------------------
    series_csv = OUTPUT_DIR / f"{BASENAME}_series.csv"
    probes_csv = OUTPUT_DIR / f"{BASENAME}_probes.csv"
    png = OUTPUT_DIR / f"{BASENAME}_series.png"
    pts = chk["pts"]
    ref = chk["ref"]
    qh_all = np.real(run["s_h"].values)
    if comm.rank == 0:
        with series_csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            head = ["n", "t_s"]
            for j in range(3):
                head += [f"probe_c{j}_K", f"ball_exact_c{j}_K", f"point_exact_c{j}_K", f"probe_qh_c{j}_K"]
            head += ["probe_c2_displaced_K"]
            for a in A_SWEEP:
                head += [f"probe_c0_a{a * 1e3:.0f}mm_K", f"ball_exact_c0_a{a * 1e3:.0f}mm_K"]
            w.writerow(head)
            for n in range(N_MAIN + 1):
                row = [n, f"{n * DT:.6f}"]
                for j in range(3):
                    row += [f"{np.real(s_cf.values[n, j]):.9e}", f"{ref[n, j]:.9e}",
                            f"{pts[n, j]:.9e}", f"{qh_all[n, j]:.9e}"]
                row += [f"{np.real(s_cf.values[n, 3]):.9e}"]
                for k, a in enumerate(A_SWEEP):
                    row += [f"{Px[n, k + 1]:.9e}", f"{ex[a][n]:.9e}"]
                w.writerow(row)
        with probes_csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["label", "cx_m", "cy_m", "cz_m", "a_m", "V_h_over_V_minus_1"])
            names = ["c0", "c1", "c2", "c2_displaced"]
            for j, c in enumerate(CENTRES + [DISPLACED]):
                w.writerow([names[j], *(f"{v:.9e}" for v in c), f"{A_PROBE:.3e}", f"{vol[j]:.9e}"])
            for k, a in enumerate(A_SWEEP):
                w.writerow([f"c0_a{a * 1e3:.0f}mm", 0.0, 0.0, 0.0, f"{a:.3e}", f"{sv[k + 1]:.9e}"])

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8))
        cols = ["tab:blue", "tab:green", "tab:red"]
        names = ["c0 (0)", "c1 (3.5 mm)", "c2 (6.5 mm)"]
        for j in range(3):
            ax1.plot(ts, ref[:, j], "-", color=cols[j], label=f"{names[j]} exact ball mean")
            ax1.plot(ts, pts[:, j], "--", color=cols[j], label=f"{names[j]} point value")
            ax1.plot(ts, np.real(s_cf.values[:, j]), "o", ms=3, color=cols[j])
        ax1.plot(ts, np.real(s_cf.values[:, 3]), "x", ms=4, color="k", label="c2 displaced 1 mm (C1)")
        ax1.set_xlabel("t [s]")
        ax1.set_ylabel("theta [K]")
        ax1.set_title("a = 2 mm probes (markers: FEM)")
        ax1.legend(fontsize=7)
        aa = np.linspace(0.5e-3, 8.0e-3, 30)
        curve = [abs(_ball_mean(0.0, a, T_STAR, q) / pt0[N_MAIN] - 1.0) for a in aa]
        ax2.loglog(aa * 1e3, curve, "-", label="exact ball / point - 1")
        ax2.loglog(np.array([a * 1e3 for a in labels]), [abs(off_fem[a]) for a in labels], "o",
                   label="FEM probe / point - 1")
        ax2.loglog(aa * 1e3, curve[len(curve) // 2] * (aa / aa[len(curve) // 2]) ** 2, ":",
                   color="gray", label="a^2 slope")
        ax2.axhline(P_BAND, color="r", ls="--", label="1 % band")
        ax2.set_xlabel("a [mm] at c0")
        ax2.set_ylabel("|P / theta(0, t*) - 1|")
        ax2.set_title("volume offset at t*")
        ax2.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(png, dpi=110)
        plt.close(fig)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf_path} -- 41 steps: P1 theta_h, theta_exact, theta_error; DG0 probe masks")
    _say(comm, f"[output] {series_csv}")
    _say(comm, f"[output] {probes_csv}")
    _say(comm, f"[output] {png}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size} "
               f"(fixture {run['t_total']:.1f} s)")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
