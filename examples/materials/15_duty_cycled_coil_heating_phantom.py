"""Example (`EX-99`): a duty-cycled coil drive heats the gel phantom -- the on/off history.

`mat:9` drives `ANS-2`'s 10 MHz F-small gel phantom with a source that never
switches off; `mat:11` switches a *uniform* source on `MAT-9`'s ball.  This example
switches the **coil-driven** source (single-leg drive k = 0, G = 1 W / P_1) on and
off -- 50 % duty, period P = t*/5, every switch on a grid point -- and shows the
time history, Duhamel's sum and the half-power stand-in.

No new physics, and nothing restated.  ``_run_duty(return_fields=True)`` is
`MAT-17` step 3's own fixture body lifted to module level (rule (a) lift, EX-99;
``_transient_s`` gained an additive ``snapshot_steps=`` keyword, default the old
literal); the bands ``L_BAND`` / ``D_BAND`` / ``E_BAND`` / ``D_CONTROL_FLOOR``
are imported by name.  The example adds two thermal runs that keep every full-step
snapshot (n = 0 ... 40), for the duty schedule and for chi = 1, so it can write
<theta>(t_n), a ParaView time series and the fields at t*.

Asserted, each printed beside ``20261006T201547Z_MAT-17-step3-W1.log:1956-1976``:
(S) the schedule, (Delta) Duhamel superposition <= ``L_BAND``, (D') chi = 1 is P1
<= ``D_BAND``, (E) every run and (E') at N and n = 36 <= ``E_BAND``, (Hist) the
history runs equal the gate's runs <= ``D_BAND``; negative control (Delta-ctl) one
Duhamel term dropped, > ``L_BAND``.  Printed, never asserted: control (d) at both
placements, (c-sched) the schedule ignored, the records at t*, the history.

Scope: 1 W normalisation (no absolute temperature), degree 1, one rung, no
perfusion; not a pulse sequence, SAR averaging, the tuned coil, 64 MHz or a
periodic steady state.

Run through the harness (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:15 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_15_duty_cycled_coil_heating_phantom_combined.xdmf``,
``..._history.xdmf`` and ``..._mean_history.csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import (
    adopt_host_ownership,
    write_xdmf_time_series,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation import test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
from tests.validation import test_mat15_coil_driven_heating_phantom as m15  # noqa: E402
from tests.validation import test_mat17_duty_cycle_coil_driven as m17s3  # noqa: E402
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import PHANTOM_CELL_TAG  # noqa: E402
from tests.validation.test_mat17_duty_cycle_coil_driven import (  # noqa: E402
    D_BAND,
    D_CONTROL_FLOOR,
    DT,
    DUTY,
    E_BAND,
    L_BAND,
    N_MAIN,
    N_ON_LAST,
    SNAPS,
    T_INF,
    T_STAR,
    _allsum,
    _chi_duty,
    _run_duty,
    _transient_s,
)
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_15_duty_cycled_coil_heating_phantom"
W1 = "20261006T201547Z_MAT-17-step3-W1.log"


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
        "m17s3.E_BAND is m12.E_BAND": m17s3.E_BAND is m12.E_BAND,
        "m17s3.D_CONTROL_FLOOR is m15.D_CONTROL_FLOOR":
            m17s3.D_CONTROL_FLOOR is m15.D_CONTROL_FLOOR,
        "m17s3._run_heated is m15._run_heated": m17s3._run_heated is m15._run_heated,
        "m17s3.T_STAR is m15.T_STAR": m17s3.T_STAR is m15.T_STAR,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    if rank0:
        print("=" * 72)
        print("EX-99 -- a duty-cycled coil drive heats the gel phantom (MAT-17 step 3, re-executed)")
        print("=" * 72)
        print("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
        print(f"[records] every anchor below is beside {W1}:1956-1976", flush=True)

    # ---- the gate's own fixture, with its fields ----------------------------
    r = _run_duty(return_fields=True)
    base, r_duty, r_one = r["base"], r["r_duty"], r["r_one"]
    sub, ft, q0, dep = r["sub"], r["ft"], r["q0"], r["dep"]
    msh, cell_tags = base["msh"], base["cell_tags"]
    v_sub = base["v_sub"]
    T_p1 = base["T_star"][0]

    # EX-57 setup figure (after the fixture, which builds the mesh); no-op unless
    # FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            CONDUCTOR_CELL_TAG: "conductor",
            2: "air",
            PHANTOM_CELL_TAG: "phantom (gel)",
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="mat:15 -- a duty-cycled coil drive heats the gel phantom (MAT-17 step 3)",
        comm=comm,
    )

    # ---- the two history runs: every full-step snapshot ---------------------
    all_steps = range(N_MAIN + 1)
    h_duty = _transient_s(sub, ft, q0, _chi_duty, snapshot_steps=all_steps)
    h_one = _transient_s(sub, ft, q0, lambda t: 1.0, snapshot_steps=all_steps)

    V1 = r_duty.temperature.function_space
    work = fem.Function(V1)

    def _l2sq(arr):
        work.x.array[:] = arr
        work.x.scatter_forward()
        return float(np.real(_allsum(comm, ufl.inner(work, work) * ufl.dx(domain=sub))))

    def _mean_theta(T):
        return float(np.real(_allsum(comm, (T - T_INF) * ufl.dx(domain=sub)))) / v_sub

    def _rel(a, b):
        return float(np.sqrt(_l2sq(a.x.array - b.x.array) / _l2sq(b.x.array - T_INF)))

    # ---- (S) the schedule ----------------------------------------------------
    chis = r_duty.source_multipliers
    chi_str = "".join(str(int(c)) if c in (0.0, 1.0) else "?" for c in chis)
    on_time = sum(chis[2:]) * DT + 0.5 * DT * sum(chis[:2])
    sched_equal = list(h_duty.source_multipliers) == list(chis)

    # ---- (Hist) the history runs are the gate's runs ------------------------
    hist_duty = _rel(h_duty.temperature, r_duty.temperature)
    hist_one = _rel(h_one.temperature, T_p1)
    hist_snaps = {n: _rel(h_duty.snapshots[n], r_duty.snapshots[n]) for n in SNAPS}

    # ---- records at t* (gate format and !r) ---------------------------------
    mean_duty = _mean_theta(r_duty.temperature)
    mean_p1 = _mean_theta(T_p1)
    ratio = mean_duty / (DUTY * mean_p1)
    n_lo, n_hi = SNAPS[0], N_ON_LAST
    m_lo = _mean_theta(r_duty.snapshots[n_lo])
    m_hi = _mean_theta(r_duty.snapshots[n_hi])
    # (c-sched) the schedule ignored: theta[1] against theta_duty at t*
    csched = float(np.sqrt(_l2sq(r_one.temperature.x.array - r_duty.temperature.x.array)
                           / _l2sq(r_duty.temperature.x.array - T_INF)))

    # ---- the mean history ----------------------------------------------------
    rows = []
    for n in all_steps:
        t_n = n * DT
        if n == 0:
            chi_txt = ""
        elif n == 1:
            chi_txt = f"{h_duty.source_multipliers[0]:g}/{h_duty.source_multipliers[1]:g}"
        else:
            chi_txt = f"{h_duty.source_multipliers[n]:g}"
        md = _mean_theta(h_duty.snapshots[n])
        mh = DUTY * _mean_theta(h_one.snapshots[n])
        rows.append((n, t_n, chi_txt, md, mh, (md / mh) if n > 0 else float("nan")))

    if rank0:
        print(f"\n[anchor] (S) chi per (half-)step ({len(chis)}) = {chi_str}  "
              f"(W1: 11111000011110000111100001111000011110000); on-time {on_time:.4f} s "
              f"(W1: 627.0000); source_multipliers of the history run equal the gate's: "
              f"{sched_equal}")
        print(f"[anchor] (Delta) Duhamel superposition {r['Delta']:.3e} <= {L_BAND:.0e}  "
              f"(W1: 3.421e-14)")
        print(f"[anchor] (D') chi=1 vs step 1 {r['D']:.3e} <= {D_BAND:.0e}  (W1: 2.119e-15)")
        worst = max(r["E"], key=r["E"].get)
        print(f"[anchor] (E) worst over {len(r['E'])} runs {r['E'][worst]:.3e} on {worst} "
              f"<= {E_BAND:.0e}  (W1: 3.743e-14); history runs {h_duty.max_balance_residual:.3e} "
              f"(duty), {h_one.max_balance_residual:.3e} (chi=1)")
        print(f"[anchor] (E') step N={N_MAIN} {r['Ep_N']:.3e} <= {E_BAND:.0e}  (W1: 3.317e-15); "
              f"step n={N_ON_LAST} {r['Ep_on']:.3e} <= {E_BAND:.0e}  (W1: 1.399e-14)")
        print(f"[anchor] (Hist) history duty run vs gate's {hist_duty:.3e}; chi=1 vs "
              f"T_star[0] {hist_one:.3e}; snapshots "
              + ", ".join(f"n={n}: {v:.3e}" for n, v in hist_snaps.items())
              + f"  (all <= {D_BAND:.0e}; predicted <= 1e-15)")
        print(f"[control] (Delta-ctl) one Duhamel term (m=9) dropped: rel L2 {r['Delta_ctl']:.6f} "
              f"> {L_BAND:.0e}  (W1: 0.298210)")
        print(f"[printed] (d) doubled source, (E') at n={N_ON_LAST} (on, chi=1) "
              f"{r['ctl_d']:.6f} (W1: 1.000000, gate asserts > {D_CONTROL_FLOOR:.0e}); at "
              f"N={N_MAIN} (off, chi=0, blind by linearity) {r['ctl_d_N']:.3e} (W1: 2.942e-15)")
        print(f"[printed] (c-sched) schedule ignored, ||theta[1] - theta_duty||/||theta_duty|| at "
              f"t* = {csched:.6f} (predicted 1.0-1.2, bracket [1.067, 3.067])")
        print(f"[record] <theta>(t*): duty {mean_duty:.9f} K vs D*P1 = {DUTY} * {mean_p1:.9f} K = "
              f"{DUTY * mean_p1:.9f} K ratio {ratio:.6f}  (W1: 0.568835514 / 0.576656825 / "
              f"0.986437)")
        print(f"[record!r] duty {mean_duty!r} P1 {mean_p1!r} ratio {ratio!r} "
              f"gate mean_ratio {r['mean_ratio']!r}")
        print(f"[record] ripple, last on/off pair: n={n_lo} {m_lo:.9f} K -> n={n_hi} (on end) "
              f"{m_hi:.9f} K -> n={N_MAIN} (off end) {mean_duty:.9f} K  (W1: 0.466250330 -> "
              f"0.584438780 -> 0.568835514)")
        print(f"[record!r] m_lo {m_lo!r} m_hi {m_hi!r}")
        print(f"[record] max|Im theta| = {r['max_imag']:.3e}  (W1: 0)")
        print(f"[record] 10 g ball theta(c_0) at t*: duty {r['ball10_duty']:.9e} K  P1 "
              f"{r['ball10_p1']:.9e} K  duty/(D*P1) = "
              f"{r['ball10_duty'] / (DUTY * r['ball10_p1']):.6f}  (W1 :1975)")
        print(f"[record!r] max_imag {r['max_imag']!r} Delta {r['Delta']!r} "
              f"Delta_ctl {r['Delta_ctl']!r} csched {csched!r}")
        print(f"[integers] sub-cell count {sub.topology.index_map(sub.topology.dim).size_global}; "
              f"chi string {chi_str}; chi list {[float(c) for c in chis]}")
        print("\n[history] n  t_n(s)  chi  <theta>_duty(K)  D*<theta>_P1(K)  ratio")
        for n, t_n, chi_txt, md, mh, rt in rows:
            print(f"[history] {n:2d} {t_n:8.2f} {chi_txt:>4s} {md:.9f} {mh:.9f} "
                  f"{rt:.6f}")
        print(f"[history!r] t* row {rows[-1][3]!r} {rows[-1][4]!r}", flush=True)

    # ---- asserts (every rank) -------------------------------------------------
    assert r["chi_binary"], "(S) a midpoint chi is not 0/1"
    assert sched_equal, "(S) history-run schedule differs from the gate's"
    assert r["Delta"] <= L_BAND, f"(Delta) {r['Delta']:.3e} > {L_BAND:.0e}"
    assert r["D"] <= D_BAND, f"(D') {r['D']:.3e} > {D_BAND:.0e}"
    for name, val in r["E"].items():
        assert val <= E_BAND, f"(E) run {name}: {val:.3e}"
    assert h_duty.max_balance_residual <= E_BAND, "(E) history duty run"
    assert h_one.max_balance_residual <= E_BAND, "(E) history chi=1 run"
    assert r["Ep_N"] <= E_BAND, f"(E') step N: {r['Ep_N']:.3e}"
    assert r["Ep_on"] <= E_BAND, f"(E') last on step: {r['Ep_on']:.3e}"
    assert hist_duty <= D_BAND, f"(Hist) duty run {hist_duty:.3e}"
    assert hist_one <= D_BAND, f"(Hist) chi=1 run {hist_one:.3e}"
    for n, v in hist_snaps.items():
        assert v <= D_BAND, f"(Hist) snapshot n={n}: {v:.3e}"
    assert r["Delta_ctl"] > L_BAND, f"(Delta-ctl) did not miss: {r['Delta_ctl']:.3e}"

    # ---- outputs -------------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)

    def _field(name, arr):
        f = fem.Function(V1, name=name)
        f.x.array[:] = arr
        f.x.scatter_forward()
        return f

    th_duty = _field("theta_duty", r_duty.temperature.x.array - T_INF)
    th_on = _field("theta_on_end", r_duty.snapshots[N_ON_LAST].x.array - T_INF)
    th_half = _field("theta_half_power", DUTY * (T_p1.x.array - T_INF))
    q0.name = "q_h_k0"  # the XDMF dataset name is the Function's, not the dict key
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None,
        {"theta_duty": th_duty, "theta_on_end": th_on, "theta_half_power": th_half,
         "q_h_k0": q0},
        comm=comm,
    )
    steps = [(n * DT, {"theta": _field("theta", h_duty.snapshots[n].x.array - T_INF)})
             for n in all_steps]
    hist_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_history", sub, None, steps, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_mean_history.csv"
        with open(csv_path, "w") as fh:
            fh.write("n,t_s,chi,mean_theta_duty_K,D_mean_theta_P1_K,ratio\n")
            for n, t_n, chi_txt, md, mh, rt in rows:
                fh.write(f"{n},{t_n!r},{chi_txt},{md!r},{mh!r},{'' if n == 0 else repr(rt)}\n")
        adopt_host_ownership(OUTPUT_DIR, comm=MPI.COMM_SELF)
        print(f"\n[output] {xdmf_path.name}: P1 theta_duty, theta_on_end (n={N_ON_LAST}), "
              f"theta_half_power (K, t* = {T_STAR:.1f} s) and DG0 q_h_k0 (W/m^3 at G = 1 W/P_1)")
        print(f"[output] {hist_path.name}: P1 theta at {len(steps)} steps, t_n = n*{DT:.4f} s")
        print(f"[output] {csv_path.name}: {len(rows)} rows")
        print("EX-99 anchors and control: PASS")


if __name__ == "__main__":
    main()
