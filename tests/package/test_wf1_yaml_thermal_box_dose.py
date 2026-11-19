"""`WF-1` step 30 (§9 item 878, §10 U139): the ``thermal_box`` YAML gains a
``thermal.dose`` readout (`MAT-20`'s ``cem43`` and six-minute trailing window)
on U125's walled-box files, and it reproduces `MAT-21` P68's dose gate.

Brief: docs/planning/chunks/WF-1.md "§9 item 878".  Real build, ``-n 2``.  The
parse node builds no mesh; the solve node runs only with ``FEM_EM_WF1_S30=1``.

Anchors (bands imported, none set): (Q-wf) ``cem43`` on the FEM point history
against ``_q_check_g``'s quadrature at ``DOSE_BAR`` — P68's gate
(``20261105T224651Z_MAT-21-P68-W1.log:78–83``, ≤ 5.374e-16), on P68's mesh
(119 808 cells, ``:58``) and Δt (70.5375 s, ``:53``); (Win-wf) the trailing mean
at ``AVG_BAR``; (I-h) / (I-d) the file's history and dose against a direct R1
in the same window at ``IDENTITY_BAND``; (R) U125's θ(P5)
(``20261105T210638Z_WF-1-s28-W1-r2.log:91``) within 5e-13; (E) ``E_BAND``.
Control (rule (e), asserted at the gate's band): (c-init-d) R0's ``cem43`` at
P5 misses R1's by > ``DOSE_BAR`` (P68's form; its D = 0.25 miss read
4.6009e-01, ``…P68-W1.log:85``).

Not claimed: P68's D = 0.25 dose values (a wall beside a duty is refused), a
dose on the ball / EM sphere / unwalled box / under a duty, a dose band, limit
or threshold, compliance, a cited PMMA ρc.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    DoseReadout,
    DutyCycle,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    UniformHeatSource,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study import TransverseQuadraticHeatSource
from fem_em_solver.study_yaml import _load, _study_from_mapping
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.dose import cem43, running_average
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_thermal_box_wall import (
    A_BOX,
    B_BOX,
    C_BOX,
    CELLS,
    D_WALL,
    DT,
    E_BAND,
    G_GRAD,
    H_GEL,
    IDENTITY_BAND,
    N_MAIN,
    R0_FILE,
    R1_FILE,
    RHO_C,
    ROBIN,
    RUNG,
    T_INF,
    K,
    Q,
    _field,
    _materials,
    _py_study,
    _rho_c_field,
    _robin_facets,
    _walled_box,
    _with,
    _without,
)
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR, T_GRID
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import _dose
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_fem import (
    GRID_BAR,
    HIST_LABELS,
    HIST_PTS,
)
from tests.validation.test_mat21_wall_pre_equilibration_fem import _read

S30_ENV = "FEM_EM_WF1_S30"
R1_DOSE = WF1 / "thermal_box_wall_gradient_m4_dose.yaml"
R0_DOSE = WF1 / "thermal_box_wall_m4_dose.yaml"
MALFORMED = WF1 / "malformed_thermal_box_dose_wall_without_thickness.yaml"
STEADY = WF1 / "thermal_box_steady_m4.yaml"
DUTY_BOX = WF1 / "thermal_box_duty_m6.yaml"
QUAD_BOX = WF1 / "thermal_box_quadratic_cn_m4.yaml"
BALL = WF1 / "thermal_ball_duty_a.yaml"
DOSE_LINES = 3  # the appended ``thermal.dose`` block: ``dose:``, ``points_m``, ``window_s``
#: U125's θ(P5) on R1 (`20261105T210638Z_WF-1-s28-W1-r2.log:91`); bound = half the last digit.
U125_THETA_P5 = 1.224895450933e00
R_BOUND = 5e-13


def _knob() -> bool:
    raw = os.environ.get(S30_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S30_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _raises(exc, fn):
    try:
        fn()
    except exc as e:
        return str(e)
    return None


def _readout() -> DoseReadout:
    return DoseReadout(tuple(map(tuple, HIST_PTS)), W_B2)


def _rel(a, b):
    return abs(a - b) / abs(b)


# ------------------------------------------------------------------ helpers for item 888


def _direct_r1_history(comm=MPI.COMM_WORLD):
    """P63's direct R1 (``test_wf1_yaml_thermal_box_wall.py:330–338``) with every full-step
    snapshot, each read by ``_read`` on ``HIST_PTS`` (one batch), then cleared.

    Returns (history T = ``T_INF`` + θ, (n_t × n_pts), times, result, cells, ok, solve s, read s).
    """
    tq = time.perf_counter()
    msh, _ns, _nw = _walled_box(RUNG, (False, False, True))
    chi, kf, _inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731
    kw = dict(
        k=kf,
        q=Q * chi,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN,
        theta=0.5,
        startup_half_steps=2,
        dt=DT,
    )
    spec_r1 = dict(
        rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad, snapshot_steps=tuple(range(N_MAIN + 1))
    )
    direct = solve_transient_conduction(msh, ft, **kw, **spec_r1)
    ts = time.perf_counter()
    hist = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
    ok_all = [True] * HIST_PTS.shape[0]
    for n in range(N_MAIN + 1):
        th, ok = _read(direct.snapshots[n], HIST_PTS)
        hist[n] = T_INF + th
        ok_all = [a and b for a, b in zip(ok_all, ok)]
    direct.snapshots.clear()
    times = DT * np.arange(N_MAIN + 1)
    cells = int(msh.topology.index_map(3).size_global)
    return hist, times, direct, cells, ok_all, ts - tq, time.perf_counter() - ts


def _dose_gates(label, T_hist, times):
    """(Q-wf) and (Win-wf) on one history: {name: bool} and the printed rows."""
    gates, rows = {}, []
    for j, lab in enumerate(HIST_LABELS):
        dose, quad, rel, _swap, miss = _q_check_g(T_hist[:, j], times)
        mx, t_mx, lo, hi, ref, wrel = _win_check_g(T_hist[:, j], times, window=W_B2)
        gates[f"Q_wf_{label}_{lab}"] = rel <= DOSE_BAR
        gates[f"Win_wf_{label}_{lab}"] = (wrel <= AVG_BAR) and (lo <= mx <= hi)
        rows.append(
            f"  (Q-wf) {label} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
            f"(DOSE_BAR {DOSE_BAR:g}); R swap miss {miss:.4e} (printed); (Win-wf) peak "
            f"{mx:.12e} C at t {t_mx:.4f} s, quad {ref:.12e}, rel {wrel:.3e} (AVG_BAR "
            f"{AVG_BAR:g}), [lo, hi] [{lo:.9e}, {hi:.9e}]"
        )
    return gates, rows


# ------------------------------------------------------------------ no mesh


def test_wf1_s30_parse_and_refusals(monkeypatch, capsys):
    """W0: custody, the spec, the parsed readout, the refusals; no mesh, no solve."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    monkeypatch.setattr(ThermalBoxGeometry, "build", stub)
    monkeypatch.setattr(ThermalBallGeometry, "build", stub)

    # (K) custody: the dose files are U125's text plus the dose block
    for f, base in ((R1_DOSE, R1_FILE), (R0_DOSE, R0_FILE)):
        lines = f.read_text().splitlines(keepends=True)
        head, tail = "".join(lines[:-DOSE_LINES]), lines[-DOSE_LINES:]
        same = head == base.read_text()
        _say(
            capsys,
            f"\n[WF-1 s30] (K) custody {f.name}: equals {base.name} without the dose "
            f"lines {same}; dose lines {[ln.rstrip() for ln in tail]}",
        )
        assert same and tail[0] == "  dose:\n", tail
        diff = _spec_diff(load_study(f), load_study(base))
        _say(capsys, f"  (K) {f.name}: spec diff vs {base.name} {diff}")
        assert diff == [], diff
        _, kind, outputs, th = _load(f)
        _, _, outputs_b, th_b = _load(base)
        assert kind == "thermal_box" and set(outputs) == {"dose"} and outputs_b == {}
        assert repr(th) == repr(th_b), (th, th_b)
        want = _readout()
        got = outputs["dose"]
        hex_ok = all(
            float(a).hex() == float(b).hex()
            for pa, pb in zip(got.points_m, want.points_m)
            for a, b in zip(pa, pb)
        )
        _say(
            capsys,
            f"  (K) {f.name}: readout {got!r}; want {want!r}; coords by float.hex "
            f"{hex_ok}; window {float(got.window_s).hex()} vs W_B2 {float(W_B2).hex()}",
        )
        assert repr(got) == repr(want) and hex_ok and len(got.points_m) == len(want.points_m)

    for f in (R1_DOSE, R0_DOSE):
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        capsys.readouterr()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}")
        assert rc_all == (0, 0)

    # (Ref) the malformed file
    f_mal = _field(lambda: load_study(MALFORMED))
    d = LOGS / "thermal_box_dose_malformed"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(
        capsys,
        f"  (Ref) {MALFORMED.name}: field {f_mal!r}; run rc {run_all}, check rc "
        f"{check_all}, out exists {exists}",
    )
    assert f_mal == "geometry.wall.thickness_m"
    assert run_all == (2, 2) and check_all == (2, 2) and not exists

    # (Ref) the variant table
    m = _mapping(R1_DOSE)
    dose_sec = m["thermal"]["dose"]
    rows = [
        (_with(_mapping(STEADY), dose_sec, "thermal", "dose"), "thermal.dose"),
        (_with(_mapping(DUTY_BOX), dose_sec, "thermal", "dose"), "thermal.dose"),
        (_with(_mapping(QUAD_BOX), dose_sec, "thermal", "dose"), "thermal.dose"),
        (
            _with(m, [[0.0, 0.0, 0.045], [0.0, 0.0, 0.05]], "thermal", "dose", "points_m"),
            "thermal.dose.points_m[1]",
        ),
        (_with(m, {"directory": "x", "probes": []}, "outputs"), "outputs.probes"),
        (_with(_mapping(BALL), dose_sec, "thermal", "dose"), "thermal.dose"),
        (_without(m, "thermal", "dose", "points_m"), "thermal.dose.points_m"),
        (_with(m, 0.0, "thermal", "dose", "window_s"), "thermal.dose.window_s"),
    ]
    got = [_field(lambda mm=mm: _study_from_mapping(mm)) for mm, _ in rows]
    for (mm, want_f), gf in zip(rows, got):
        _say(capsys, f"  (Ref) variant -> {gf!r} (want {want_f!r})")
    assert got == [wf for _, wf in rows], got

    # (Ref) the API refusals, each before any mesh
    from tests.package.test_wf15_study_lossy_sphere import _study as _em_study

    cool = RobinCooling(K, H_GEL, T_INF)
    tr = TransientSpec(DT, N_MAIN, RHO_C)
    dr = _readout()
    half = (A_BOX, B_BOX, C_BOX)
    box = ThermalBoxGeometry(half, tuple(_mapping(STEADY)["geometry"]["divisions"]))
    api = {
        "steady box": _raises(
            ValueError, lambda: Study(box, None, UniformHeatSource(Q)).temperature(cool, dose=dr)
        ),
        "ball": _raises(
            ValueError,
            lambda: Study(ThermalBallGeometry(0.05, 0.01), None, UniformHeatSource(Q)).temperature(
                cool, transient=tr, dose=dr
            ),
        ),
        "EM sphere": _raises(ValueError, lambda: _em_study(degree=1).temperature(cool, dose=dr)),
        "beside a duty": _raises(
            ValueError,
            lambda: Study(box, None, UniformHeatSource(Q)).temperature(
                cool, transient=tr, duty=DutyCycle(564.3, 0.5), dose=dr
            ),
        ),
        "quadratic source": _raises(
            ValueError,
            lambda: Study(box, None, TransverseQuadraticHeatSource(Q)).temperature(
                cool, transient=tr, dose=dr
            ),
        ),
        "point outside the gel": _raises(
            ValueError,
            lambda: _py_study().temperature(
                cool, transient=tr, dose=DoseReadout(((0.0, 0.0, 0.05),))
            ),
        ),
        "empty points": _raises(ValueError, lambda: DoseReadout(())),
        "window 0": _raises(ValueError, lambda: DoseReadout(((0.0, 0.0, 0.0),), 0.0)),
    }
    for k_, v in api.items():
        _say(capsys, f"  (Ref) API {k_}: {v!r}")
    assert all(v is not None for v in api.values()), api
    _say(capsys, f"  build calls {calls['n']}")
    assert calls["n"] == 0


# ------------------------------------------------------------------ the solve


def test_wf1_s30_dose_reproduces_p68_gate(capsys):
    """W1: (Grid), (Cells), (Q-wf), (Win-wf), (I-h)/(I-d), (R), (E), CLI = API, (c-init-d)."""
    if not _knob():
        pytest.skip(f"{S30_ENV} unset: the walled-box dose solve is environment-selected")
    assert not IS_COMPLEX, "the walled-box records are real-build"
    comm = MPI.COMM_WORLD
    assert comm.size == 2, comm.size
    t0 = time.perf_counter()

    h_d, t_d, direct, cells_d, ok_d, s_solve, s_read = _direct_r1_history(comm)
    _say(
        capsys,
        f"\n[WF-1 s30] direct R1: cells {cells_d}; (E) "
        f"{direct.max_balance_residual:.3e}; theta(P5, t_end) {h_d[-1, HIST_LABELS.index('P5')] - T_INF:+.12e}"
        f" K; solve {s_solve:.1f} s, {N_MAIN + 1} snapshot reads {s_read:.1f} s; batch "
        f"{HIST_LABELS} {HIST_PTS.tolist()}",
    )

    runs = {}
    for key, f in (("R1", R1_DOSE), ("R0", R0_DOSE), ("R1_nodose", R1_FILE)):
        tq = time.perf_counter()
        out = LOGS / f"thermal_box_dose_{key}"
        _fresh_dir(out, comm)
        run = run_study_file(f, output_dir=out, comm=comm)
        on_disk = None
        if comm.rank == 0:
            on_disk = json.loads((out / "summary.json").read_text())
        on_disk = comm.bcast(on_disk, root=0)
        runs[key] = dict(run=run, disk=on_disk, s=time.perf_counter() - tq)
        tm = run.summary["temperature"]
        _say(
            capsys,
            f"  file {f.name}: cells {run.summary['cells']}; (E) "
            f"{tm['max_balance_residual']:.3e}; dose key {'dose' in run.summary}; "
            f"{runs[key]['s']:.1f} s (solve + reads)",
        )

    tq = time.perf_counter()
    py = _py_study().temperature(
        RobinCooling(K, H_GEL, T_INF),
        transient=TransientSpec(DT, N_MAIN, RHO_C, initial_gradient_k_per_m=G_GRAD),
        dose=_readout(),
    )
    s_py = time.perf_counter() - tq
    cells_py = int(py.mesh.topology.index_map(3).size_global)
    _say(
        capsys,
        f"  Python Study R1 + dose: cells {cells_py}; (E) "
        f"{py.max_balance_residual:.3e}; {s_py:.1f} s (solve + reads)",
    )

    dz = {k: runs[k]["run"].summary["dose"] for k in ("R1", "R0")}
    dz["py"] = py.dose
    hist = {k: np.asarray(v["temperature_c"], dtype=float) for k, v in dz.items()}
    times = {k: np.asarray(v["times_s"], dtype=float) for k, v in dz.items()}

    gates = {"valid_points": all(ok_d)}
    # (Grid)
    grid = {k: float(np.max(np.abs(t - T_GRID)) / DT) for k, t in times.items()}
    _say(
        capsys,
        f"  (Grid) max |t - T_GRID| / DT {grid} (GRID_BAR {GRID_BAR:g}); hex t_end "
        f"{ {k: float(t[-1]).hex() for k, t in times.items()} } vs direct "
        f"{float(t_d[-1]).hex()}",
    )
    for k, v in grid.items():
        gates[f"Grid_{k}"] = (times[k].shape == T_GRID.shape) and v <= GRID_BAR
    # (Cells)
    gates["Cells"] = (
        cells_d == CELLS
        and cells_py == CELLS
        and all(v["run"].summary["cells"] == CELLS for v in runs.values())
    )
    # (Q-wf) / (Win-wf)
    for k in ("R1", "R0", "py"):
        g, rows = _dose_gates(k, hist[k], times[k])
        gates.update(g)
        for r in rows:
            _say(capsys, r)
    # (I-h) / (I-d)
    d_direct = [_dose(h_d[:, j] - T_INF) for j in range(len(HIST_LABELS))]
    for k in ("R1", "py"):
        ih = float(np.max(np.abs(hist[k] - h_d) / np.abs(h_d)))
        idd = [_rel(dz[k]["cem43_s"][j], d_direct[j]) for j in range(len(HIST_LABELS))]
        _say(
            capsys,
            f"  (I-h) {k} history vs direct max rel {ih:.3e}; (I-d) cem43 rel per point "
            f"{[f'{x:.3e}' for x in idd]} (IDENTITY_BAND {IDENTITY_BAND:.0e}); D {k} "
            f"{[f'{x:.6e}' for x in dz[k]['cem43_s']]} vs direct "
            f"{[f'{x:.6e}' for x in d_direct]}",
        )
        gates[f"I_h_{k}"] = ih <= IDENTITY_BAND
        gates[f"I_d_{k}"] = all(x <= IDENTITY_BAND for x in idd)
    # (R)
    j5 = HIST_LABELS.index("P5")
    th5 = hist["R1"][-1, j5] - T_INF
    _say(
        capsys,
        f"  (R) R1 theta(P5, t_end) {th5:+.12e} K vs U125 {U125_THETA_P5:+.12e}, |d| "
        f"{abs(th5 - U125_THETA_P5):.3e} (bound {R_BOUND:g})",
    )
    gates["R"] = abs(th5 - U125_THETA_P5) <= R_BOUND
    # (E)
    gates["E_direct"] = float(direct.max_balance_residual) <= E_BAND
    gates["E_python"] = float(py.max_balance_residual) <= E_BAND
    for k, v in runs.items():
        gates[f"E_{k}"] = v["run"].summary["temperature"]["max_balance_residual"] <= E_BAND
    # CLI = API by repr
    for k in ("R1", "R0", "py"):
        z = dz[k]
        tt, TT = np.asarray(z["times_s"]), np.asarray(z["temperature_c"])
        ok_c, ok_w = [], []
        for j in range(len(HIST_LABELS)):
            ok_c.append(repr(float(cem43(tt, TT[:, j]))) == repr(z["cem43_s"][j]))
            _m, mx, t_mx = running_average(tt, TT[:, j], window=z["six_minute"]["window_s"])
            ok_w.append(
                repr(float(mx)) == repr(z["six_minute"]["peak_c"][j])
                and repr(float(t_mx)) == repr(z["six_minute"]["t_max_s"][j])
            )
        _say(
            capsys,
            f"  (CLI = API) {k}: cem43 by repr {ok_c}; six-minute peak / t_max by repr "
            f"{ok_w}; six_minute {z['six_minute']}",
        )
        gates[f"repr_{k}"] = all(ok_c) and all(ok_w)
    disk_diffs = {}
    for k_, v in runs.items():
        mine = json.loads(json.dumps(v["run"].summary))
        disk_diffs[k_] = sorted(
            kk
            for kk in set(mine) | set(v["disk"])
            if kk != "elapsed_s" and mine.get(kk) != v["disk"].get(kk)
        )
    _say(
        capsys,
        f"  summary vs summary.json differing keys (elapsed_s excluded) per rank "
        f"{comm.gather(disk_diffs, root=0)}",
    )
    gates["summary_on_disk"] = all(dd == [] for dd in disk_diffs.values())
    gates["summary_dose_keys"] = (
        "dose" in runs["R1"]["run"].summary
        and "dose" in runs["R0"]["run"].summary
        and "dose" not in runs["R1_nodose"]["run"].summary
    )
    # (c-init-d), rule (e): asserted at the gate's band
    d1, d0 = dz["R1"]["cem43_s"][j5], dz["R0"]["cem43_s"][j5]
    miss = abs(d0 / d1 - 1.0)
    _say(
        capsys,
        f"  (c-init-d) P5: R0 D {d0:.6e} s vs R1 D {d1:.6e} s, miss {miss:.4e} "
        f"(asserted > {DOSE_BAR:g}; predicted, printed only: |miss| > 0 {miss > 0.0} and at least "
        f"of order U125's theta rel miss 1.493e-01: {miss >= 1.493e-01}); signed R0/R1 - 1 "
        f"{d0 / d1 - 1.0:+.4e}",
    )
    j_a = HIST_LABELS.index("ADJ_PT")
    _say(
        capsys,
        f"  D (printed) R1 P5 {d1:.6e} ADJ_PT {dz['R1']['cem43_s'][j_a]:.6e}; R0 P5 "
        f"{d0:.6e} ADJ_PT {dz['R0']['cem43_s'][j_a]:.6e} s",
    )
    gates["c_init_d"] = miss > DOSE_BAR
    _say(capsys, f"  gates {gates}; elapsed {time.perf_counter() - t0:.1f} s")
    assert all(gates.values()), gates
