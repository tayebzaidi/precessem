"""Example (`EX-85`): a duty-cycled source on the gel ball -- what the average hides.

`mat:4`, `mat:6`, `mat:8` all heat the 10 cm gel ball with a source that is
always on. This one switches the same uniform source 0/1 (`MAT-17` step 2's
``source_schedule``): case A (period P = 1200 s, duty D = 0.5) and case B
(P = 0.2 t_R, D = 0.25), `MAT-9`'s ball (k = 0.6, Q = 2000 W/m^3, rho c = 4.18e6,
Robin h = 10 W/m^2/K, Bi = 5/3, T_inf = T(0) = 20 C), 0.005 m mesh, P1, CN +
Rannacher. Three things are shown: the on/off history against the closed form;
the periodic steady state's centre peak against the steady state of the
time-averaged power (closed form, step 1, printed); and what a continuous source
of the same average power hides, as P / tau_1 goes from 0.063 (case A) to 0.73
(case B).

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every helper
and band is imported from ``tests/validation/test_mat17_duty_cycle_fem.py``
(`m17`), ``test_mat17_duty_cycle_closed_form.py`` (`m17cf`) and the complex-safe
twins of ``test_mat13_transient_perfusion_sphere.py`` (`m13`); `MAT-9`'s own
``_interior_error`` / ``_rel_diff`` call ``float()`` on a complex value and are
never used. No gate edit, no ``src/`` edit.

Asserted (imported bands): (S0) sum chi*dt = on-time on both scheduled runs
(inside ``m17._census``); (iii) every run <= ``BALANCE_BAND``; (i) the gate's
comparison re-run in the complex build, cases A and B <= ``INTERIOR_BAND``; (H)
two new mid-phase instants (case A step 90 = 13 500 s, case B step 4 = DP/2);
(AVG) the full-power run against the duty-1 closed form (step response).
Negative controls, asserted as band misses: (c1) the always-on run against case
A's duty comparand; (c-avg B) the average-power stand-in D*(T_full - T_inf)
against case B's duty comparand.

Printed, never gated: P / tau_1; the closed-form periodic peak / trough / mean;
the hidden peak at A step 92 and B step 8; the case-A average-power stand-in.

Run it::

    ./scripts/run_examples.sh -e mat:11 -n 2 -t 420

Output (on the ball mesh): ``examples/materials/paraview_output/
materials_11_duty_cycled_gel_ball_history_combined.xdmf`` (time series),
``..._tend_combined.xdmf`` and ``..._centre_history.csv``.
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np
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

from tests.validation import test_mat13_transient_perfusion_sphere as m13  # noqa: E402
from tests.validation import test_mat17_duty_cycle_closed_form as m17cf  # noqa: E402
from tests.validation import test_mat17_duty_cycle_fem as m17  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_11_duty_cycled_gel_ball"
RES = 0.005
DT_A = 150.0
LEFT = 1.0 - 1e-12  # step 1's left-limit idiom at a switch

#: Gate records, printed beside this run and never asserted.
#: W1b = 20261003T171113Z_MAT-17-step2-W1b.log ; S1 = 20261003T010113Z_MAT-17-step1-W1.log
REC = {
    "e_a": "4.3870e-04", "e_b": "1.0683e-03", "c1": "9.7087e-01",
    "th_a": (3.053003, 3.050426), "th_b": (1.537399, 1.535862),
    "bal": 2.823e-14,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _field(V, name, values):
    f = fem.Function(V, name=name)
    f.x.array[:] = values
    f.x.scatter_forward()
    return f


def _periodic(P, D):
    """Closed-form periodic centre peak / trough / mean (step 1's idiom)."""
    peak = float(m17cf._modal(0.0, D * P * LEFT, P, D, periodic=True))
    trough = float(m17cf._modal(0.0, P * LEFT, P, D, periodic=True))
    return peak, trough, D * m17cf.TH_SS0


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # ---- structural checks (the gate modules have not drifted apart) ---------------
    assert m17.CASES is m17cf.CASES
    assert m17.INTERIOR_BAND == m13.INTERIOR_BAND and m17.BALANCE_BAND == m13.BALANCE_BAND
    assert (m17.R, m17.K, m17.Q, m17.H, m17.RHO_C, m17.T_INF, m17.ROBIN_TAG) == (
        m13.R, m13.K, m13.Q, m13.H, m13.RHO_C, m13.T_INF, m13.ROBIN_TAG)
    assert (m17.R, m17.K, m17.RHO_C, m17.Q, m17.H, m17.T_INF) == (
        m17cf.R, m17cf.K, m17cf.RHO_C, m17cf.Q, m17cf.H, m17cf.T_INF)
    BAND, BAL = m13.INTERIOR_BAND, m13.BALANCE_BAND
    T_INF = m17.T_INF
    P_A, D_A, P_B, D_B = m17.P_A, m17.D_A, m17.P_B, m17.D_B
    T_END_A, T_END_B, DT_B = m17.T_END_A, m17.T_END_B, m17.DT_B

    _say(comm, "=" * 78)
    _say(comm, "EX-85 -- a duty-cycled source on the gel ball: what the average hides")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] structural checks pass; bands imported: interior {BAND:.0e}, "
               f"balance {BAL:.0e}")

    # ---- time constants (printed) --------------------------------------------------
    tau1, lam1 = m17cf.TAU1, m17cf.LAM1
    _say(comm, f"[time] lambda_1 = {lam1:.6f}, tau_1 = {tau1:.6e} s (S1 log :141); "
               f"P/tau_1: case A {P_A / tau1:.4f}, case B {P_B / tau1:.4f} "
               f"(predicted 0.0626 / 0.727)")

    # ---- periodic steady state, closed form (printed) -------------------------------
    per = {}
    for lab, (P, D) in (("A", (P_A, D_A)), ("B", (P_B, D_B))):
        peak, trough, mean = _periodic(P, D)
        jstar = math.ceil(30.0 * tau1 / P)
        per[lab] = (peak, trough, mean)
        _say(comm, f"[periodic] case {lab} (closed form, step 1 gated): centre mean "
                   f"D*theta_ss = {mean:.6f} K, peak {peak:.6f} K, trough {trough:.6f} K, "
                   f"ripple p-p {peak - trough:.6f} K; peak/(D theta_ss) - 1 = "
                   f"{100 * (peak / mean - 1):+.3f} % (S1 log :144 A +1.174 %, :239 B "
                   f"+20.440 %); j* = ceil(30 tau_1/P) = {jstar} (S1 :142 / :237: 479 / 42)")
    rr0, tt0 = m17cf._grid(P_A, D_A)
    jc = None
    for j in range(0, 400):
        d = (m17cf._modal(rr0, tt0 + j * P_A, P_A, D_A)
             - m17cf._modal(rr0, tt0 + j * P_A, P_A, D_A, periodic=True))
        if float(np.max(np.abs(d))) <= 1e-3 * D_A * m17cf.TH_SS0:
            jc = j
            break
    _say(comm, f"[periodic] case A cycles to within 1e-3 D theta_ss of periodic: {jc} "
               f"(S1 :234: 111)")

    # ---- mesh + setup figure ------------------------------------------------------------
    msh, cell_tags, ft = m17._ball(RES)
    n_cells = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] MAT-9's _ball({RES}): {n_cells} cells (W1b: 151548)")
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel ball (Robin h = 10, Bi = 5/3; duty-cycled uniform source)"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:11 -- 10 cm Robin gel ball (MAT-9 fixture), duty-cycled source",
        comm=comm,
    )

    # ---- the four runs ----------------------------------------------------------------------
    n_a = m17._n(T_END_A, DT_A)
    n_b = m17._n(T_END_B, DT_B)
    assert n_a == 94 and n_b == 20
    snaps_a, snaps_b = tuple(range(0, n_a + 1)), tuple(range(0, n_b + 1))
    s_ad = m17._run(msh, ft, DT_A, n_a, 0.5, 2, schedule=m17._chi(P_A, D_A),
                    snapshot_steps=snaps_a)
    s_af = m17._run(msh, ft, DT_A, n_a, 0.5, 2, snapshot_steps=snaps_a)
    s_bd = m17._run(msh, ft, DT_B, n_b, 0.5, 2, schedule=m17._chi(P_B, D_B),
                    snapshot_steps=snaps_b)
    s_bf = m17._run(msh, ft, DT_B, n_b, 0.5, 2, snapshot_steps=snaps_b)
    _say(comm, f"[runs] A-duty {s_ad.wall:.1f} s, A-full {s_af.wall:.1f} s, B-duty "
               f"{s_bd.wall:.1f} s, B-full {s_bf.wall:.1f} s")

    # (S0): asserts internally
    m17._census(s_ad, DT_A, 2, m17.ON_TIME_A, "A-duty")
    m17._census(s_bd, DT_B, 2, m17.ON_TIME_B, "B-duty")
    bals = {"A-duty": s_ad.max_balance_residual, "A-full": s_af.max_balance_residual,
            "B-duty": s_bd.max_balance_residual, "B-full": s_bf.max_balance_residual}
    worst = max(bals, key=bals.get)
    _say(comm, f"[balance] (iii) worst per-step balance {bals[worst]:.3e} on '{worst}' of 4 "
               f"(record {REC['bal']:.3e}, W1b :107)  ASSERTED <= {BAL:.0e}")

    # ---- the five closed-form interpolations (9.5 s each on this mesh) ---------------------------
    ex_a = m17._exact_duty(msh, T_END_A, P_A, D_A)
    ex_b = m17._exact_duty(msh, T_END_B, P_B, D_B)
    ex_a90 = m17._exact_duty(msh, 90 * DT_A, P_A, D_A)
    ex_b4 = m17._exact_duty(msh, 4 * DT_B, P_B, D_B)
    ex_a_full = m17._exact_duty(msh, T_END_A, P_A, 1.0)

    V1 = fem.functionspace(msh, ("Lagrange", 1))

    def avg_T(T_full, D, name):
        """T_inf + D (T_full - T_inf): exact for the linear theta-scheme (MAT-17 (L))."""
        return _field(T_full.function_space, name, T_INF + D * (T_full.x.array - T_INF))

    # ---- (i): the gate's comparison re-run in the complex build ------------------------------------------
    e_a = m13._interior_error(msh, s_ad.temperature, ex_a)
    e_b = m13._interior_error(msh, s_bd.temperature, ex_b)
    th0_a, th0_b = m17._theta_centre(s_ad.temperature), m17._theta_centre(s_bd.temperature)
    cf_a = float(m17._closed(0.0, T_END_A, P_A, D_A)[0])
    cf_b = float(m17._closed(0.0, T_END_B, P_B, D_B)[0])
    _say(comm, f"\n[(i)] case A CN+Rannacher dt=150 N=94: interior rel L2 {e_a:.4e} (record "
               f"{REC['e_a']}, W1b :98; string {'equal' if f'{e_a:.4e}' == REC['e_a'] else 'DIFFERENT'}); "
               f"theta_h(0) {th0_a:.6f} K vs closed form {cf_a:.6f} K (W1b {REC['th_a'][0]} vs "
               f"{REC['th_a'][1]})  ASSERTED <= {BAND:.0e}")
    _say(comm, f"[(i)] case B CN+Rannacher dt={DT_B:.4f} N=20: interior rel L2 {e_b:.4e} (record "
               f"{REC['e_b']}, W1b :101; string {'equal' if f'{e_b:.4e}' == REC['e_b'] else 'DIFFERENT'}); "
               f"theta_h(0) {th0_b:.6f} K vs closed form {cf_b:.6f} K (W1b {REC['th_b'][0]} vs "
               f"{REC['th_b'][1]})  ASSERTED <= {BAND:.0e}")

    # ---- (H): two new mid-phase instants ---------------------------------------------------------------
    e_h_a = m13._interior_error(msh, s_ad.snapshots[90], ex_a90)
    e_h_b = m13._interior_error(msh, s_bd.snapshots[4], ex_b4)
    _say(comm, f"[(H)] case A cycle 11 mid-on (step 90, t = 13500 s): interior rel L2 {e_h_a:.4e} "
               f"(predicted [1e-4, 3e-3])  ASSERTED <= {BAND:.0e}")
    _say(comm, f"[(H)] case B cycle 0 mid-on (step 4, t = DP/2 = {4 * DT_B:.2f} s): interior rel "
               f"L2 {e_h_b:.4e} (predicted [2e-4, 5e-3])  ASSERTED <= {BAND:.0e}")

    # ---- (AVG): the full-power run is the duty-1 step response -----------------------------------------------
    e_avg = m13._interior_error(msh, s_af.temperature, ex_a_full)
    _say(comm, f"[(AVG)] A-full vs the duty-1 closed form at t_end: interior rel L2 {e_avg:.4e} "
               f"(predicted ~4e-4; MAT-9 CN read 4.6327e-4 at t* = 13933 s)  ASSERTED <= {BAND:.0e}")

    # ---- negative controls ------------------------------------------------------------------------------------------
    e_c1 = m13._interior_error(msh, s_af.temperature, ex_a)
    Tavg_b = avg_T(s_bf.temperature, D_B, "T_avg_B")
    e_cavg = m13._interior_error(msh, Tavg_b, ex_b)
    _say(comm, f"[control] (c1) always-on run vs case A's duty comparand: interior rel L2 "
               f"{e_c1:.4e} (record {REC['c1']}, W1b :104; string "
               f"{'equal' if f'{e_c1:.4e}' == REC['c1'] else 'DIFFERENT'})  ASSERTED > {BAND:.0e}")
    _say(comm, f"[control] (c-avg B) average-power stand-in vs case B's duty comparand at "
               f"t_end: interior rel L2 {e_cavg:.4e} (predicted [0.15, 0.6]); ASSERTED > "
               f"{BAND:.0e}; factor {e_cavg / BAND:.1f} x band (printed only)")

    # ---- the case-A stand-in at t_end (printed) -----------------------------------------------------------------------
    Tavg_a = avg_T(s_af.temperature, D_A, "T_avg_A")
    e_sa = m13._interior_error(msh, Tavg_a, ex_a)
    _say(comm, f"[stand-in] case A, D*(T_full - T_inf) at t_end (mid-off, cycle 11) vs the "
               f"duty comparand: interior rel L2 {e_sa:.4e} (predicted [5e-3, 4e-2]; "
               f"printed only, straddles the band)")

    # ---- the centre history and the hidden peaks --------------------------------------------------------------------------------
    t0 = time.perf_counter()
    m17._theta_centre(s_ad.snapshots[0])
    t_first = time.perf_counter() - t0
    _say(comm, f"[centre] first collective _theta_centre call {t_first:.3f} s")
    full_rows = t_first <= 1.0
    steps_a_csv = list(range(0, n_a + 1, 2))
    if not full_rows:
        keep = set()
        for j in (0, 1, 10, 11):
            keep.update(range(8 * j, 8 * j + 8, 2))
        steps_a_csv = [s for s in steps_a_csv if s in keep]
        _say(comm, "[centre] CSV restricted to cycles 0, 1, 10, 11 (call > 1 s)")
    rows = []
    for case, P, D, dt, sd, sf, steps in (
        ("A", P_A, D_A, DT_A, s_ad, s_af, steps_a_csv),
        ("B", P_B, D_B, DT_B, s_bd, s_bf, list(range(0, n_b + 1))),
    ):
        for n in steps:
            t = n * dt
            h = m17._theta_centre(sd.snapshots[n])
            f = m17._theta_centre(sf.snapshots[n])
            cf = 0.0 if n == 0 else float(m17._closed(0.0, t * LEFT, P, D)[0])
            rows.append((case, n, t, h, cf, D * f, h - D * f))
    if comm.rank == 0:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        with (OUTPUT_DIR / f"{BASENAME}_centre_history.csv").open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["case", "step", "t_s", "theta_h_0_K", "closed_form_left_limit_K",
                         "D_theta_full_0_K", "duty_minus_average_K"])
            for r in rows:
                wr.writerow([r[0], r[1], f"{r[2]:.4f}", f"{r[3]:.9e}", f"{r[4]:.9e}",
                             f"{r[5]:.9e}", f"{r[6]:.9e}"])

    peak_rows = {}
    for case, n in (("A", 92), ("B", 8)):
        P, D, dt, sd, sf = ((P_A, D_A, DT_A, s_ad, s_af) if case == "A"
                            else (P_B, D_B, DT_B, s_bd, s_bf))
        t = n * dt
        h = m17._theta_centre(sd.snapshots[n])
        f = m17._theta_centre(sf.snapshots[n])
        cf = float(m17._closed(0.0, t * LEFT, P, D)[0])
        peak_rows[case] = (h, cf, D * f, h / (D * f) - 1.0)
        _say(comm, f"[hidden] case {case} step {n} (t = {t:.1f} s, end of the on phase): "
                   f"theta_h(0) {h:.6f} K, closed form (left limit) {cf:.6f} K, D*theta_full(0) "
                   f"{D * f:.6f} K; peak/avg - 1 = {100 * (h / (D * f) - 1):+.3f} % "
                   f"(predicted {'[3, 7] %' if case == 'A' else '300 % = 1/D - 1'})")
    d_b8 = m13._rel_diff(msh, s_bd.snapshots[8], s_bf.snapshots[8])
    _say(comm, f"[hidden] case B first on-phase is a full-power exposure: ||theta_duty - "
               f"theta_full||/||theta_full|| at step 8 = {d_b8:.3e} (predicted <= 1e-14)")

    # ---- assertions ---------------------------------------------------------------------------------------------------------------------
    assert bals[worst] <= BAL, f"(iii) {worst} {bals[worst]:.3e}"
    assert e_a <= BAND, f"(i) A {e_a:.4e}"
    assert e_b <= BAND, f"(i) B {e_b:.4e}"
    assert e_h_a <= BAND, f"(H) A {e_h_a:.4e}"
    assert e_h_b <= BAND, f"(H) B {e_h_b:.4e}"
    assert e_avg <= BAND, f"(AVG) {e_avg:.4e}"
    assert e_c1 > BAND, f"(c1) {e_c1:.4e} did not miss"
    assert e_cavg > BAND, f"(c-avg B) {e_cavg:.4e} did not miss"

    # ---- outputs (CG1 exports are interpolants, nothing asserted) -----------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frames = []
    for n in range(0, n_a + 1, 2):
        d = s_ad.snapshots[n].x.array - T_INF
        a = D_A * (s_af.snapshots[n].x.array - T_INF)
        frames.append((n * DT_A, {
            "theta_duty": _field(V1, "theta_duty", d),
            "theta_avg_power": _field(V1, "theta_avg_power", a),
            "theta_hidden": _field(V1, "theta_hidden", d - a),
        }))
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_history_combined", msh, None, frames, comm=comm,
    )

    fields = {}
    for lab, s_d, s_f, D, ex in (("A", s_ad, s_af, D_A, ex_a), ("B", s_bd, s_bf, D_B, ex_b)):
        th_h = _field(V1, f"theta_h_{lab}", s_d.temperature.x.array - T_INF)
        th_s = fem.Function(V1, name=f"theta_series_{lab}")
        th_s.interpolate(ex)
        th_s.x.scatter_forward()
        fields[f"theta_h_{lab}"] = th_h
        fields[f"theta_series_{lab}"] = th_s
        fields[f"theta_err_{lab}"] = _field(V1, f"theta_err_{lab}", th_h.x.array - th_s.x.array)
        fields[f"theta_avg_{lab}"] = _field(V1, f"theta_avg_{lab}",
                                            D * (s_f.temperature.x.array - T_INF))
    chi_int = m13._interior_indicator(msh)
    mask = fem.Function(chi_int.function_space, name="interior_mask")
    mask.x.array[:] = chi_int.x.array
    mask.x.scatter_forward()
    fields["interior_mask"] = mask
    tend_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_tend_combined", msh, None, fields, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {series_path} -- P1 theta_duty, theta_avg_power, theta_hidden at "
               f"case-A steps 0, 2, ..., 94 ({len(frames)} frames)")
    _say(comm, f"[output] {tend_path} -- t_end fields for cases A and B + interior_mask")
    _say(comm, f"[output] {OUTPUT_DIR / (BASENAME + '_centre_history.csv')} ({len(rows)} rows)")

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
