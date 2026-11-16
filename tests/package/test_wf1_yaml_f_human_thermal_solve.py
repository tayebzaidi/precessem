"""`WF-1` step 10b (§9 item 431b): `MAT-15` P17's ball mean and `MAT-20`'s six-minute maximum through the YAML.

``run_study_file`` on the tracked tuned-copper F-human thermal file (item
431a) at ``-n 8`` solves once and writes the transient Robin-cooled series
of the phantom; the ccw 10 g ball mean at c*_θ, the volume mean and the
360 s running maximum are read against item 380's series JSON (brief:
docs/planning/chunks/WF-1.md, "§9 item 431", part 431b).  Every band and
record is imported; none is set here.  A counting call-through wrapper on
``fem_em_solver.study._study_solve_tuned_quadrature`` reads 1 after every
read; ``build_tuned_quadrature`` is never called.

Anchors (asserted):
  (A0)    summary cells == ``HOLE_CELL_RECORD`` within ``CELL_COUNT_BAND``;
          one solve; ``run.result is run.study.result``.
  (S_rep) the summary's ``s_node`` against the node JSON <= ``REUSE_REPRODUCTION_RTOL``.
  (V)     the phantom submesh volume against ``V_PHANTOM_H9P_RECORD_M3`` <= ``V_BAND``;
          n_sub == the JSON's.
  (B)     ∫q against G·P_phantom,ccw (``mean_sar`` on the run's field) <= ``B_BAND``.
  (E)     the thermal balance residual <= ``E_BAND``.
  (T*)    θ̄(c*_θ, t*) against the JSON's last ball sample <= ``REUSE_REPRODUCTION_RTOL``.
  (Ser)   times to 1e-12 s; ball and volume-mean series elementwise within
          ``REUSE_REPRODUCTION_RTOL`` of the series maximum; the bundle JSON
          equals the summary on rank 0.
  (6min)  the six-minute peak against ``running_average`` on the JSON series
          <= ``REUSE_REPRODUCTION_RTOL``, t_max exact; summary peak against an
          in-test ``running_average`` on its own series <= ``AVG_BAR``.

Negative controls (rule (e)): (N1) the feed-bisector image M(c*_θ) misses by
> ``REUSE_REPRODUCTION_RTOL`` (asserted, backed: P17 13.21 %); (N2) the window
one sample short misses by > the band (asserted, backed: MAT-20 1.36 %);
(N3) the doubled source misses by > ``D_CONTROL_FLOOR`` (asserted; the
factor 1.000 is predicted, printed).

Opt-in: ``FEM_EM_WF1_S10B=1`` (unset ⇒ the solve nodes skip; any other value
raises); width 8 asserted when set.  Scope: degree-1 EM, 64 MHz, κ off, one
copper hole fixture, one gel phantom, ccw only.  Not claimed: the console
script, the cw sense, CEM43, an absolute temperature or compliance figure,
degree 2, perfusion.
"""

from __future__ import annotations

import dataclasses
import json
import os
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import study as study_mod
from fem_em_solver.post.sar import mean_sar
from fem_em_solver.study_yaml import _load, run_study_file
from fem_em_solver.thermal.dose import running_average

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.package.test_wf1_study_yaml import LOGS
from tests.package.test_wf1_yaml_f_human_thermal import FHT_FILE
from tests.package.test_wf16_study_tuned_f_human import _cx
from tests.package.test_wf16_study_tuned_f_human_b1_sar import RECORD_RANK_WIDTH
from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import D_CONTROL_FLOOR, V_BAND
from tests.validation.test_mat15_f_human_tuned_heating import (
    SERIES_JSON,
    V_PHANTOM_H9P_RECORD_M3,
)
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_wf10_f_human_tuned_b1_sar import NODE_JSON
from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

S10B_ENV = "FEM_EM_WF1_S10B"
REPO = Path(__file__).resolve().parents[2]
H12C_NPZ = "logs/wf12/h12c-cu-s32-a2.npz"
LOGGED_THETA_TSTAR_K = 4.926110850e-02  # 20261011T183420Z_MAT-20-step3-W1.log (printed only)
LOGGED_SIX_MIN_K = 4.277005e-02  # 20261011T183717Z_MAT-20-step3-W2.log (printed only)
LOGGED_VOL_SIX_MIN_K = 3.070e-02  # MAT-20 step 3 volume-mean six-minute peak (printed only)


def _knob() -> bool:
    raw = os.environ.get(S10B_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S10B_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-1 s10b] {msg}", flush=True)


def _inputs():
    """The three gitignored inputs, read / checked on rank 0 and broadcast."""
    comm = MPI.COMM_WORLD
    out = None
    if comm.rank == 0:
        missing = [p for p in (SERIES_JSON, NODE_JSON, H12C_NPZ) if not (REPO / p).is_file()]
        series = node = None
        if not missing:
            series = json.loads((REPO / SERIES_JSON).read_text())
            node = json.loads((REPO / NODE_JSON).read_text())
            for p in (SERIES_JSON, NODE_JSON, H12C_NPZ):
                st = os.stat(REPO / p)
                print(f"[WF-1 s10b] input {p}: {st.st_size} bytes, mtime "
                      f"{time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(st.st_mtime))}",
                      flush=True)
        out = (missing, series, node)
    return comm.bcast(out, root=0)


def test_wf1_fht_solve_knob_and_inputs():
    """No mesh: the knob, the inputs present (absent ⇒ FAIL), the file's G and c*_θ."""
    _knob()  # a bad value raises here too
    assert FHT_FILE.exists(), FHT_FILE
    missing, series, _node = _inputs()
    assert not missing, f"inputs absent: {missing} — stop (BLOCKED, rule (d))"
    _study, kind, _outputs, fh = _load(FHT_FILE)
    _p0(f"kind {kind}; source_scale {fh.cooling.source_scale.hex()} vs G "
        f"{float(series['G_per_w']).hex()}; centres {fh.balls['centres']} vs c*_theta "
        f"{series['c_star_m']}; json ranks {series['ranks']}, n_sub {series['n_sub']}, cells "
        f"{series['cells']}, {len(series['times_s'])} samples")
    assert fh.cooling.source_scale == float(series["G_per_w"])
    assert [list(map(float, c)) for c in fh.balls["centres"]] == [list(series["c_star_m"])]
    assert series["ranks"] == RECORD_RANK_WIDTH
    assert series["cells"] == HOLE_CELL_RECORD


@pytest.fixture(scope="module")
def solved():
    if not _knob():
        pytest.skip(f"{S10B_ENV} unset: the YAML F-human thermal solve is a heavy -n 8 window, "
                    "environment-selected")
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, f"records are -n {RECORD_RANK_WIDTH}; got {comm.size}"
    missing, series, node = _inputs()  # read and broadcast before the run
    assert not missing, f"inputs absent: {missing} — stop (BLOCKED, rule (d))"
    _study, _kind, _outputs, fh = _load(FHT_FILE)

    calls = {"n": 0}
    real = study_mod._study_solve_tuned_quadrature

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(study_mod, "_study_solve_tuned_quadrature", counting)
        comm.Barrier()
        t0 = time.perf_counter()
        run = run_study_file(FHT_FILE, output_dir=LOGS / "f_human_tuned_thermal", comm=comm)
        comm.Barrier()
        elapsed = time.perf_counter() - t0
        res = run.result
        _p0(f"run_study_file: {elapsed:.1f} s at -n {comm.size} (build "
            f"{res.times.get('build', 0.0):.1f} s, node pair {res.times.get('node_pair', 0.0):.1f} s); "
            f"solve calls {calls['n']}")
        yield {"run": run, "series": series, "node": node, "fh": fh, "calls": calls,
               "comm": comm, "elapsed": elapsed}


def _red(comm, form):
    return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))


@complex_only
def test_wf1_fht_cells_once_and_source(solved):
    """(A0), (S_rep), (V), (B), (E); η printed."""
    run, comm, series, node = solved["run"], solved["comm"], solved["series"], solved["node"]
    res, summary, th = run.result, run.summary, run.thermal
    cell_rel = summary["cells"] / HOLE_CELL_RECORD - 1.0
    s = np.array([[_cx(v) for v in row] for row in summary["s_node"]], dtype=np.complex128)
    s_json = np.array([[_cx(v) for v in row] for row in node["s_node"]], dtype=np.complex128)
    r_s = _residual(s, s_json)
    # (V): the phantom submesh of the run's mesh, the same helper the branch uses (collective).
    sub, _ids = study_mod._tag_submesh(res.mesh, res.cell_tags, run.study.geometry.phantom_tag)
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = _red(comm, fem.form(one * ufl.dx(domain=sub)))
    r_v = abs(v_sub / V_PHANTOM_H9P_RECORD_M3 - 1.0)
    # (B): P_phantom,ccw by mean_sar on the run's ccw field (collective reduction inside).
    f = res.fields["ccw"]
    tagged = mean_sar(f.e_complex, sigma=f.sigma_field, rho=PHANTOM_RHO_KG_PER_M3,
                      cell_tags=res.cell_tags, subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
    p_ph = float(np.real(tagged["dissipated_power_w"]))
    G = float(series["G_per_w"])
    tp = summary["temperature"]
    r_b = abs(tp["source_power_w"] - G * p_ph) / (G * p_ph)
    eta = p_ph * G  # = P_phantom / P_acc,ccw (G = 1 / P_acc,ccw)
    _p0(f"(A0) cells {summary['cells']} vs record {HOLE_CELL_RECORD} (rel {cell_rel:.3e}, band "
        f"{CELL_COUNT_BAND}, ASSERTED); thermal cells {th.cells}; solve calls "
        f"{solved['calls']['n']} (== 1 ASSERTED)")
    _p0(f"(S_rep) summary s_node vs node JSON residual {r_s:.6e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED; prior 4.287849e-12)")
    _p0(f"(V) V_sub {v_sub:.12e} m^3 vs record {V_PHANTOM_H9P_RECORD_M3:.12e}: rel {r_v:.3e} (band "
        f"{V_BAND:g}, ASSERTED; prior 7.749e-14); n_sub {th.n_sub} vs JSON {series['n_sub']}")
    _p0(f"(B) source_power_w {tp['source_power_w']:.15e} W vs G*P_phantom,ccw {G * p_ph:.15e} W "
        f"(P_phantom {p_ph:.15e} W, JSON {series['p_phantom_ccw_w']:.15e}): rel {r_b:.3e} (band "
        f"{B_BAND:g}, ASSERTED; prior 0)")
    _p0(f"(E) max balance residual {tp['max_balance_residual']:.3e} (band {E_BAND:g}, ASSERTED; "
        f"P17 <= 3.3e-10)")
    _p0(f"eta = P_phantom / P_acc,ccw {eta:.6f} (printed; P17 0.2821)")
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert summary["cells"] == int(res.cells) == th.cells
    assert run.result is run.study.result
    assert r_s <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_s:.6e}"
    assert r_v <= V_BAND, f"(V) {r_v:.3e}"
    assert th.n_sub == series["n_sub"], (th.n_sub, series["n_sub"])
    assert r_b <= B_BAND, f"(B) {r_b:.3e}"
    assert tp["max_balance_residual"] <= E_BAND, tp["max_balance_residual"]
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]


@complex_only
def test_wf1_fht_ball_reproduces_p17(solved):
    """(T*): θ̄(c*_θ, t*) against the JSON; the ball volume printed."""
    run, series = solved["run"], solved["series"]
    tp = run.summary["temperature"]
    assert len(tp["balls"]) == 1
    ball = tp["balls"][0]
    th_star = ball["theta_k"][-1]
    th_json = float(series["theta_ball10g_k"][-1])
    r_t = abs(th_star / th_json - 1.0)
    _p0(f"(T*) theta-bar(c*_theta, t* = {tp['times_s'][-1]:.1f} s) {th_star:.12e} K vs JSON "
        f"{th_json:.12e} K: rel {r_t:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED; predicted "
        f"<= 1e-10); logged {LOGGED_THETA_TSTAR_K:.9e} K")
    _p0(f"ball radius {ball['radius_m']:.15e} m vs JSON a10 {series['a10_m']:.15e}; volume "
        f"{ball['volume_m3']:.15e} m^3 vs JSON {series['ball_volume_m3']:.15e} (rel "
        f"{abs(ball['volume_m3'] / series['ball_volume_m3'] - 1.0):.3e}, printed); sense "
        f"{tp['sense']}; centre {ball['centre']}")
    assert tp["sense"] == "ccw"
    assert r_t <= REUSE_REPRODUCTION_RTOL, f"(T*) {r_t:.3e}"
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]


@complex_only
def test_wf1_fht_series_and_six_minute(solved):
    """(Ser), (6min); the volume-mean six-minute peak printed."""
    run, series = solved["run"], solved["series"]
    tp = run.summary["temperature"]
    t_sum = np.array(tp["times_s"])
    t_json = np.array(series["times_s"])
    assert t_sum.shape == t_json.shape, (t_sum.shape, t_json.shape)
    d_t = float(np.max(np.abs(t_sum - t_json)))
    b_sum = np.array(tp["balls"][0]["theta_k"])
    b_json = np.array(series["theta_ball10g_k"])
    m_sum = np.array(tp["volume_mean_k"])
    m_json = np.array(series["theta_mean_k"])
    r_ball = float(np.max(np.abs(b_sum - b_json)) / np.max(np.abs(b_json)))
    r_mean = float(np.max(np.abs(m_sum - m_json)) / np.max(np.abs(m_json)))
    _m, peak_json, tmax_json = running_average(t_json, b_json, window=W_B2)
    six = tp["six_minute"]
    peak, tmax = six["balls"][0]["peak_k"], six["balls"][0]["t_max_s"]
    r_6 = abs(peak / float(peak_json) - 1.0)
    _m2, peak_own, _t_own = running_average(t_sum, b_sum, window=W_B2)
    r_own = abs(peak - float(peak_own)) / abs(float(peak_own))
    _m3, peak_vol_json, _tv = running_average(t_json, m_json, window=W_B2)
    _p0(f"(Ser) {t_sum.size} samples; max|t - t_json| {d_t:.3e} s (<= 1e-12 ASSERTED); ball "
        f"series max|d|/max {r_ball:.3e}, volume-mean {r_mean:.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED); sample 0 ball {b_sum[0]!r} mean {m_sum[0]!r}")
    _p0(f"(6min) window {six['window_s']} s: peak {peak:.12e} K at {tmax} s vs JSON series "
        f"{float(peak_json):.12e} K at {float(tmax_json)} s: rel {r_6:.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED; t_max exact ASSERTED); logged "
        f"{LOGGED_SIX_MIN_K:.6e} K at 1254 s")
    _p0(f"(6min) API vs in-test running_average on the summary's own series rel {r_own:.3e} (band "
        f"{AVG_BAR:g}, ASSERTED)")
    _p0(f"volume-mean six-minute peak {six['volume_mean']['peak_k']:.9e} K at "
        f"{six['volume_mean']['t_max_s']} s (JSON series {float(peak_vol_json):.9e}; MAT-20 "
        f"{LOGGED_VOL_SIX_MIN_K:.3e} K; printed only)")
    assert six["window_s"] == W_B2
    assert d_t <= 1e-12, d_t
    assert r_ball <= REUSE_REPRODUCTION_RTOL, f"(Ser) ball {r_ball:.3e}"
    assert r_mean <= REUSE_REPRODUCTION_RTOL, f"(Ser) mean {r_mean:.3e}"
    assert r_6 <= REUSE_REPRODUCTION_RTOL, f"(6min) {r_6:.3e}"
    assert tmax == float(tmax_json), (tmax, tmax_json)
    assert r_own <= AVG_BAR, f"(6min) API vs operator {r_own:.3e}"
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]


@complex_only
def test_wf1_fht_controls(solved):
    """(N1) mirror ball, (N2) short window, (N3) doubled source — no re-solve."""
    run, comm, series, fh = solved["run"], solved["comm"], solved["series"], solved["fh"]
    tp = run.summary["temperature"]
    th_json = float(series["theta_ball10g_k"][-1])
    c_star = np.array([series["c_star_m"]], dtype=float)
    m_c = [float(x) for x in _mirror_xy(c_star, float(run.result.az_mirror))[0]]
    t0 = time.perf_counter()
    n1 = run.study.temperature(fh.cooling, transient=fh.transient,
                               balls={"mass_g": fh.balls["mass_g"], "centres": [m_c]},
                               window_s=fh.window_s, sense="ccw", comm=comm)
    t_n1 = time.perf_counter() - t0
    r_n1 = abs(n1.balls[0]["theta_k"][-1] / th_json - 1.0)
    dt = float(series["dt_s"])
    _m, peak_json, _t = running_average(series["times_s"], series["theta_ball10g_k"], window=W_B2)
    _m2, peak_short, _t2 = running_average(tp["times_s"], tp["balls"][0]["theta_k"],
                                           window=W_B2 - dt)
    r_n2 = abs(float(peak_short) / float(peak_json) - 1.0)
    doubled = dataclasses.replace(fh.cooling, source_scale=2.0 * fh.cooling.source_scale)
    t0 = time.perf_counter()
    n3 = run.study.temperature(doubled, transient=fh.transient, balls=fh.balls,
                               window_s=fh.window_s, sense="ccw", comm=comm)
    t_n3 = time.perf_counter() - t0
    r_n3 = abs(n3.balls[0]["theta_k"][-1] / th_json - 1.0)
    _p0(f"(N1) M(c*_theta) = {m_c} (az_mirror {run.result.az_mirror:.12f} deg): theta-bar(t*) "
        f"{n1.balls[0]['theta_k'][-1]:.9e} K, miss {r_n1:.6e} = {r_n1 / REUSE_REPRODUCTION_RTOL:.3e}x "
        f"band (> band ASSERTED; P17 mis-paired 13.21 %); {t_n1:.1f} s")
    _p0(f"(N2) window {W_B2 - dt} s peak {float(peak_short):.9e} K vs {float(peak_json):.9e} K: "
        f"miss {r_n2:.6e} = {r_n2 / REUSE_REPRODUCTION_RTOL:.3e}x band (> band ASSERTED; MAT-20 "
        f"(Sh) 1.36 %)")
    _p0(f"(N3) source_scale 2G: theta-bar(t*) {n3.balls[0]['theta_k'][-1]:.9e} K, miss "
        f"{r_n3:.6f} (> {D_CONTROL_FLOOR:g} ASSERTED; predicted 1.000); {t_n3:.1f} s")
    assert r_n1 > REUSE_REPRODUCTION_RTOL, f"(N1) blind comparator {r_n1:.3e}"
    assert r_n2 > REUSE_REPRODUCTION_RTOL, f"(N2) blind comparator {r_n2:.3e}"
    assert r_n3 > D_CONTROL_FLOOR, f"(N3) blind comparator {r_n3:.3e}"
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]


@complex_only
def test_wf1_fht_summary_json(solved):
    """(Ser) the written bundle equals the summary on rank 0; no Touchstone."""
    run, comm = solved["run"], solved["comm"]
    summary = run.summary
    with open(run.paths["summary"]) as fh:
        on_disk = json.load(fh)
    differing = sorted(k for k in set(on_disk) | set(summary) if on_disk.get(k) != summary.get(k))
    if comm.rank != 0:
        print(f"[WF-1 s10b] (J) rank {comm.rank} differing keys {differing}", flush=True)
    _p0(f"(J) rank 0 differing keys {differing} (== [] ASSERTED); family {summary['family']}; keys "
        f"{sorted(summary)}; temperature keys {sorted(summary['temperature'])}; paths "
        f"{sorted(run.paths)}; window run_study_file {solved['elapsed']:.1f} s")
    if comm.rank == 0:
        assert on_disk == summary, differing
    assert summary["family"] == "f_human_birdcage"
    assert "touchstone" not in summary and "sar" not in summary and "b1_plus" not in summary
    assert summary["temperature"]["source_scale"] == float(solved["series"]["G_per_w"])
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]
