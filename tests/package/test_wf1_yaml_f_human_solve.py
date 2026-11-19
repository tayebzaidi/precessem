"""`WF-1` step 7b (§9 item 390b, chain U43 part 2): ``run_study_file`` on the
tuned copper F-human YAML file at ``-n 8`` reproduces items 382 / 382b's
records from its own summary.

The tracked file (item 390a) parses to item 382's ``Study`` bit for bit; here
``run_study_file`` solves it once and the summary's node S, active ``Y_in``,
10 g SAR lattice and ``|B₁⁺|`` points are read against the stored node JSON
and the registered H9′ records.  Every band and record is imported; none is
set here.  The solve is never re-entered: a counting call-through wrapper on
``fem_em_solver.study._study_solve_tuned_quadrature`` reads 1 after every
test-side read.  ``build_tuned_quadrature`` is never called.

Anchors (asserted):
  (A0)    summary cells == ``HOLE_CELL_RECORD`` within ``CELL_COUNT_BAND`` and
          == ``run.result.cells``; one solve; ``run.result is run.study.result``.
  (S_rep) the summary's ``s_node`` against the node JSON <= ``REUSE_REPRODUCTION_RTOL``.
  (Y)     per sense, the summary's ``y_in`` against the JSON's, same band.
  (S10)   the 18 ccw lattice values, SAR10g(c*) and cw at M(c*) against the
          records, same band; the lattice argmax is ``H9P_CSTAR_LABEL``.
  (iii-a) the mirror at c* <= ``C4_COVARIANCE_BAND``.
  (B)     mean ``|B₁⁺|``_ccw over the 41 sample points on the run's Study
          against ``H9P_MEAN_B1_PLUS_CCW_RECORD_T``, ``REUSE_REPRODUCTION_RTOL``;
          (Sense) > 1; (M_B) <= ``C4_COVARIANCE_BAND``; (B2′) summary == API
          within ``EXACT_IDENTITY_RTOL``.
  (J)     the written JSON equals the summary on rank 0; no Touchstone.

Negative controls (rule (e), asserted, backed): (N1) the mis-paired sense >
``C4_COVARIANCE_BAND`` (prior 9.5476e-01); (N2) the summary's ``s_node``
against the H17 node JSON (118.2000 pF) misses by > ``REUSE_REPRODUCTION_RTOL``
(prior 4.266875e-02).  Factors printed, not asserted.

Printed only: (Δc), the transcribed centres against the run-derived ones.

Opt-in: ``FEM_EM_WF1_S7B=1`` (unset ⇒ the solve nodes skip; any other value
raises); width 8 asserted when set.  Scope: degree 1, 64 MHz, κ off, one
copper hole fixture, one phantom.  Not claimed: the console script, (iii-b),
the 18 cw(M) rows, an absolute ``|B₁⁺|`` or SAR, a global SAR maximum,
compliance, ``f₀`` / Q / tuning, degree 2.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import study as study_mod
from fem_em_solver.study_yaml import _load, run_study_file
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.package.test_wf1_study_yaml import LOGS
from tests.package.test_wf1_yaml_f_human import FH_FILE, _expected_centres
from tests.package.test_wf16_study_tuned_f_human import H17_NODE_JSON, _c_f, _cx
from tests.package.test_wf16_study_tuned_f_human_b1_sar import (
    H9P_CSTAR_LABEL,
    H9P_MEAN_B1_PLUS_CCW_RECORD_T,
    H9P_N_SAMPLE_POINTS,
    H9P_SAR10G_CSTAR_RECORD_W_PER_KG,
    H9P_SAR10G_TABLE_RECORD,
    RECORD_RANK_WIDTH,
)
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    _relative_l2,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy, _read_senses
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_wf7_f_human_sar10g_hotspot import _phantom_geometry
from tests.validation.test_wf10_f_human_tuned_b1_sar import NODE_JSON
from tests.validation.test_wf10_f_human_tuned_sar10g import AZ17_RECORD_DEG, _lattice
from tests.validation.test_wf12_f_human_copper_matrix import COPPER_SIGMA, HOLE_CELL_RECORD

S7B_ENV = "FEM_EM_WF1_S7B"
N_LATTICE = 18  # the file's centres 0..17 are _lattice's rows, 18 is M(c*)


def _knob() -> bool:
    raw = os.environ.get(S7B_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S7B_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-1 s7b] {msg}", flush=True)


def test_wf1_fh_solve_knob_and_imports():
    """No mesh: the knob, the file's lattice against the records' labels."""
    _knob()  # a bad value raises here too
    assert FH_FILE.exists(), FH_FILE
    _study, _kind, outputs, _cooling = _load(FH_FILE)
    centres = [tuple(c) for c in outputs["sar"]["centres"]]
    expected = _expected_centres()
    assert len(centres) == N_LATTICE + 1 and centres == expected
    labels = [row[0] for row in _lattice(AZ17_RECORD_DEG, 0.0)]
    assert labels == [row[0] for row in H9P_SAR10G_TABLE_RECORD]
    i_star = labels.index(H9P_CSTAR_LABEL)
    twin = [i for i in range(N_LATTICE) if centres[i] == centres[N_LATTICE]]
    _p0(
        f"(K-c) 19 centres parse as transcribed; c* = row {i_star} ({H9P_CSTAR_LABEL}); "
        f"M(c*) bitwise equal to lattice row(s) {twin} "
        f"({[labels[i] for i in twin]}): the 56.25 deg mirror maps az17 onto az17+90"
    )
    assert twin == [labels.index("r=45 dz=+0 az17+90")], twin
    rec = {row[0]: row[1:] for row in H9P_SAR10G_TABLE_RECORD}
    assert rec[H9P_CSTAR_LABEL][0] == H9P_SAR10G_CSTAR_RECORD_W_PER_KG


@pytest.fixture(scope="module")
def solved():
    if not _knob():
        pytest.skip(
            f"{S7B_ENV} unset: the YAML F-human solve is a heavy -n 8 window, "
            "environment-selected"
        )
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, f"records are -n {RECORD_RANK_WIDTH}; got {comm.size}"
    # Copied (provenance: tests/package/test_wf16_study_tuned_f_human.py:192–209,
    # item 382's ``tuned`` fixture; the npz cell read is not needed here): the
    # node JSON read on rank 0 and broadcast before any solve.
    rec, h17 = None, None
    if comm.rank == 0:
        assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: stop (BLOCKED, rule (d))"
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
        if os.path.exists(H17_NODE_JSON):
            with open(H17_NODE_JSON) as fh:
                h17 = json.load(fh)
        st = os.stat(NODE_JSON)
        print(
            f"[WF-1 s7b] node JSON {NODE_JSON}: {st.st_size} bytes, mtime "
            f"{time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(st.st_mtime))}; "
            f"c_f {rec['c_f']!r}",
            flush=True,
        )
    rec, h17 = comm.bcast((rec, h17), root=0)
    assert rec["c_f"] == _c_f(), (rec["c_f"], _c_f())

    calls = {"n": 0}
    real = study_mod._study_solve_tuned_quadrature

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(study_mod, "_study_solve_tuned_quadrature", counting)
        comm.Barrier()
        t0 = time.perf_counter()
        run = run_study_file(FH_FILE, output_dir=LOGS / "f_human_tuned", comm=comm)
        comm.Barrier()
        elapsed = time.perf_counter() - t0
        res = run.result
        _p0(
            f"run_study_file: {elapsed:.1f} s at -n {comm.size} (build "
            f"{res.times.get('build', 0.0):.1f} s, node pair {res.times.get('node_pair', 0.0):.1f} s); "
            f"solve calls {calls['n']}"
        )
        yield {"run": run, "rec": rec, "h17": h17, "calls": calls, "comm": comm, "elapsed": elapsed}


def _summary_s(summary):
    return np.array([[_cx(v) for v in row] for row in summary["s_node"]], dtype=np.complex128)


@complex_only
def test_wf1_fh_yaml_cells_centres_once(solved):
    """(A0); (Δc) printed."""
    run, comm = solved["run"], solved["comm"]
    res, summary = run.result, run.summary
    cell_rel = summary["cells"] / HOLE_CELL_RECORD - 1.0
    centre, _r, _hh, _lo, _hi = _phantom_geometry(res.mesh, res.cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    lattice = _lattice(float(res.az_feeds[0]), z_c)
    derived = [tuple(float(v) for v in row[4]) for row in lattice]
    i_star = [row[0] for row in lattice].index(H9P_CSTAR_LABEL)
    derived.append(
        tuple(float(v) for v in _mirror_xy(np.array([derived[i_star]]), float(res.az_mirror))[0])
    )
    file_c = np.array(summary["sar"]["centres"], dtype=float)
    d_c = float(np.max(np.abs(file_c - np.array(derived))))
    _p0(
        f"(A0) cells {summary['cells']} vs record {HOLE_CELL_RECORD} (rel {cell_rel:.3e}, band "
        f"{CELL_COUNT_BAND}, ASSERTED), result.cells {res.cells}; solve calls {solved['calls']['n']} "
        f"(== 1 ASSERTED); feeds az {res.az_feeds[0]:.12f} / {res.az_feeds[1]:.12f}, az_mirror "
        f"{res.az_mirror:.12f}"
    )
    _p0(
        f"(Δc) max |file centre - run-derived centre| {d_c:.3e} m over 19 (printed only; "
        f"z_c {z_c:.6e} m)"
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert summary["cells"] == int(res.cells)
    assert run.result is run.study.result
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]


@complex_only
def test_wf1_fh_yaml_node_reproduces_h9p(solved):
    """(S_rep), (Y), (N2)."""
    run, rec = solved["run"], solved["rec"]
    summary = run.summary
    s = _summary_s(summary)
    s_json = np.array([[_cx(v) for v in row] for row in rec["s_node"]], dtype=np.complex128)
    r_s = _residual(s, s_json)
    _p0(
        f"(S_rep) summary s_node vs node JSON residual {r_s:.6e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED; prior 4.265317e-12)"
    )
    for i in range(2):
        for j in range(2):
            _p0(
                f"(S_rep) S[{i},{j}] summary {s[i, j]:.12e} json {s_json[i, j]:.12e} "
                f"|d| {abs(s[i, j] - s_json[i, j]):.3e}"
            )
    r_y = {}
    for sense in ("ccw", "cw"):
        y = np.array([_cx(p) for p in summary["y_in"][sense]])
        y_json = np.array([_cx(p) for p in rec["y_in"][sense]])
        r_y[sense] = float(np.linalg.norm(y - y_json) / np.linalg.norm(y_json))
        _p0(
            f"(Y) [{sense}] summary {y[0]:.10e} / {y[1]:.10e} S; JSON {y_json[0]:.10e} / "
            f"{y_json[1]:.10e}; rel {r_y[sense]:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED)"
        )
    h17 = solved["h17"]
    r_n2 = None
    if h17 is not None and "s_node" in h17:
        s17 = np.array([[_cx(v) for v in row] for row in h17["s_node"]], dtype=np.complex128)
        r_n2 = _residual(s, s17)
        _p0(
            f"(N2) summary s_node vs H17 node JSON at {h17.get('c_f')!r} F residual {r_n2:.6e} = "
            f"{r_n2 / REUSE_REPRODUCTION_RTOL:.3e}x band (> band ASSERTED; prior 4.266875e-02)"
        )
    else:
        _p0("(N2) not run: H17 node JSON absent or without s_node; (N1) carries rule (e)")
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]
    assert r_s <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_s:.6e}"
    for sense, r in r_y.items():
        assert r <= REUSE_REPRODUCTION_RTOL, f"(Y) {sense} {r:.3e}"
    if r_n2 is not None:
        assert r_n2 > REUSE_REPRODUCTION_RTOL, f"(N2) blind comparator {r_n2:.6e}"


@complex_only
def test_wf1_fh_yaml_sar_reproduces_h9p(solved):
    """(S10), (iii-a); the z = z_c twelve-point mirror printed."""
    run = solved["run"]
    sar = run.summary["sar"]
    ccw = sar["values"]["ccw"]
    cw = sar["values"]["cw"]
    assert len(ccw) == N_LATTICE + 1 and len(cw) == N_LATTICE + 1
    labels = [row[0] for row in _lattice(AZ17_RECORD_DEG, 0.0)]
    rec = {row[0]: row[1:] for row in H9P_SAR10G_TABLE_RECORD}
    sar_c = np.array([d["averaged_sar_w_per_kg"] for d in ccw])
    sar_w = np.array([d["averaged_sar_w_per_kg"] for d in cw])
    r_ccw = np.array([abs(sar_c[i] / rec[lab][0] - 1.0) for i, lab in enumerate(labels)])
    i_star = int(np.argmax(sar_c[:N_LATTICE]))
    key_star = labels[i_star]
    sar_star = float(sar_c[i_star])
    r_star = abs(sar_star - H9P_SAR10G_CSTAR_RECORD_W_PER_KG) / H9P_SAR10G_CSTAR_RECORD_W_PER_KG
    cw_m = float(sar_w[N_LATTICE])
    cw_m_rec = rec[H9P_CSTAR_LABEL][1]
    r_cw_m = abs(cw_m - cw_m_rec) / cw_m_rec
    mirror_star = abs(cw_m - sar_star) / sar_star
    _p0(f"radius_m {sar['radius_m']:.9e}, mass_g {sar['mass_g']}")
    for i, lab in enumerate(labels):
        _p0(
            f"  [{i:2d}] {lab:<22s} ccw {sar_c[i]:.9e} (rec {rec[lab][0]:.9e}, rel {r_ccw[i]:.2e}) "
            f"cw {sar_w[i]:.9e} (printed) mass {ccw[i]['mass_kg'] * 1e3:.6f} g"
            f"{'   <- c*' if i == i_star else ''}"
        )
    _p0(
        f"  [18] M(c*)                 cw {cw_m:.9e} (rec cw(M c*) {cw_m_rec:.9e}, rel "
        f"{r_cw_m:.3e}); ccw {sar_c[N_LATTICE]:.9e} (printed); mass "
        f"{cw[N_LATTICE]['mass_kg'] * 1e3:.6f} g"
    )
    twin = [i for i in range(N_LATTICE) if sar["centres"][i] == sar["centres"][N_LATTICE]]
    for i in twin:
        _p0(
            f"  M(c*) is bitwise lattice row [{i}] {labels[i]}: cw[18] - cw[{i}] = "
            f"{cw_m - sar_w[i]:.3e}, ccw[18] - ccw[{i}] = {sar_c[N_LATTICE] - sar_c[i]:.3e} "
            f"(printed)"
        )
    # z = z_c row: az17+d <-> az17+90-d pairs +0<->+90, +180<->+270 at each radius.
    idx = {lab: i for i, lab in enumerate(labels)}
    pairs = {"+0": "+90", "+90": "+0", "+180": "+270", "+270": "+180"}
    row12 = [lab for lab in labels if "dz=+0 " in lab]
    ref = np.array([sar_c[idx[lab]] for lab in row12])
    img = np.array(
        [sar_w[idx[lab.rsplit("+", 1)[0] + pairs["+" + lab.rsplit("+", 1)[1]]]] for lab in row12]
    )
    l2_12 = _relative_l2(img, ref, np.ones(ref.size, dtype=bool))
    _p0(f"z = z_c twelve-point mirror cw(Mx) vs ccw(x) rel L2 {l2_12 * 100:.4f}% (printed only)")
    _p0(
        f"(S10) ccw max rel vs records {r_ccw.max():.3e} (band {REUSE_REPRODUCTION_RTOL:g}, "
        f"ASSERTED); c* = {key_star} (record {H9P_CSTAR_LABEL}); SAR10g(c*) {sar_star:.9e} W/kg, "
        f"rel {r_star:.3e}; cw(M c*) rel {r_cw_m:.3e}"
    )
    _p0(
        f"(iii-a) mirror at c* {mirror_star * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%, "
        f"ASSERTED; record 0.6413%)"
    )
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]
    assert r_ccw.max() <= REUSE_REPRODUCTION_RTOL, f"(S10) ccw {r_ccw.max():.3e}"
    assert key_star == H9P_CSTAR_LABEL, key_star
    assert r_star <= REUSE_REPRODUCTION_RTOL, f"(S10) c* {r_star:.3e}"
    assert r_cw_m <= REUSE_REPRODUCTION_RTOL, f"(S10) cw(M c*) {r_cw_m:.3e}"
    assert mirror_star <= C4_COVARIANCE_BAND, f"(iii-a) {mirror_star:.4e} — never re-band"


@complex_only
def test_wf1_fh_yaml_b1_reproduces_h9p(solved):
    """(B), (Sense), (M_B), (B2′); (N1) asserted."""
    run, comm = solved["run"], solved["comm"]
    res, st = run.result, run.study
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})
    az_m = float(res.az_mirror)
    vals, valid = st.b1_plus_at(points, drives=("ccw",), comm=comm)["ccw"]
    vals, valid = np.asarray(vals), np.asarray(valid, dtype=bool)
    b1 = st.b1_plus(comm=comm)
    ccw_x = _read_senses(b1["ccw"], points)
    cw_mx = _read_senses(b1["cw"], _mirror_xy(points, az_m))
    mask = np.logical_and.reduce([valid, ccw_x[2], cw_mx[2]])
    n_pts, n_valid, n_mask = int(points.shape[0]), int(np.count_nonzero(valid)), int(mask.sum())
    plus_x, minus_x = ccw_x[0], ccw_x[1]
    mean_b = float(np.mean(vals[valid]))
    r_b = abs(mean_b - H9P_MEAN_B1_PLUS_CCW_RECORD_T) / H9P_MEAN_B1_PLUS_CCW_RECORD_T
    sense_ratio = float(np.mean(plus_x[mask]) / np.mean(minus_x[mask]))
    mirror = _relative_l2(cw_mx[1], plus_x, mask)
    mispaired = _relative_l2(cw_mx[0], plus_x, mask)
    # (B2′): the summary's two points against the API on the same Study.
    bp = run.summary["b1_plus"]
    pts2 = np.array(bp["points"], dtype=float)
    api = st.b1_plus_at(pts2, comm=comm)
    d2, scale, ok = 0.0, 0.0, True
    for sense, (v_api, ok_api) in api.items():
        v_api = np.asarray(v_api, dtype=float)
        v_sum = np.asarray(bp["values"][sense], dtype=float)
        d2 = max(d2, float(np.max(np.abs(v_sum - v_api))))
        scale = max(scale, float(np.max(np.abs(v_api))))
        ok = ok and all(bp["valid"][sense]) and bool(np.all(ok_api))
    r_b2 = d2 / scale
    _p0(
        f"(B) {n_pts} sample points, {n_valid} valid, {n_mask} valid on every image; mean|B1+|_ccw "
        f"{mean_b:.9e} T vs record {H9P_MEAN_B1_PLUS_CCW_RECORD_T:.6e} T: rel {r_b:.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED)"
    )
    _p0(f"(Sense) mean|B1+|_ccw / mean|B1-|_ccw {sense_ratio:.4f} (> 1 ASSERTED; record 21.5155)")
    _p0(
        f"(M_B) feed-bisector mirror (az {az_m:.6f} deg) rel L2 {mirror:.4e} (band "
        f"{C4_COVARIANCE_BAND:g}, ASSERTED; record 1.9647e-03)"
    )
    _p0(
        f"(N1) mis-paired |B1+|_cw(Mx) vs |B1+|_ccw(x) {mispaired:.4e} = "
        f"{mispaired / C4_COVARIANCE_BAND:.1f}x band (> band ASSERTED; prior 9.5476e-01)"
    )
    _p0(
        f"(B2') summary b1_plus at {pts2.shape[0]} points vs b1_plus_at: max|d| / max|B| "
        f"{r_b2:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED); valid all {ok}; values "
        f"{ {k: [f'{x:.9e}' for x in v] for k, v in bp['values'].items()} }"
    )
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]
    assert n_pts == H9P_N_SAMPLE_POINTS and n_valid == n_pts, (n_pts, n_valid)
    assert r_b <= REUSE_REPRODUCTION_RTOL, f"(B) {r_b:.3e}"
    assert sense_ratio > 1.0, f"(Sense) {sense_ratio:.4f}"
    assert mirror <= C4_COVARIANCE_BAND, f"(M_B) {mirror:.4e} — never re-band"
    assert mispaired > C4_COVARIANCE_BAND, f"(N1) blind comparator {mispaired:.4e}"
    assert r_b2 <= EXACT_IDENTITY_RTOL, f"(B2') {r_b2:.3e}"
    assert ok, "(B2') a point is invalid"


@complex_only
def test_wf1_fh_yaml_summary_json(solved):
    """(J) the written bundle."""
    run, comm = solved["run"], solved["comm"]
    summary = run.summary
    with open(run.paths["summary"]) as fh:
        on_disk = json.load(fh)
    differing = sorted(k for k in set(on_disk) | set(summary) if on_disk.get(k) != summary.get(k))
    if comm.rank != 0:
        print(f"[WF-1 s7b] (J) rank {comm.rank} differing keys {differing}", flush=True)
    _p0(
        f"(J) rank 0 differing keys {differing} (== [] ASSERTED); family {summary['family']}; "
        f"c_in_f {summary['excitation']['c_in_f'].hex()} vs {_c_f().hex()}; conductor "
        f"{summary['conductor']}; keys {sorted(summary)}; paths {sorted(run.paths)}; window "
        f"run_study_file {solved['elapsed']:.1f} s"
    )
    if comm.rank == 0:
        assert on_disk == summary, differing
    assert "touchstone" not in summary and "touchstone" not in run.paths
    assert "z_matrix" not in summary
    assert summary["family"] == "f_human_birdcage"
    assert summary["excitation"]["c_in_f"].hex() == _c_f().hex()
    assert summary["conductor"]["sigma_s_per_m"] == COPPER_SIGMA == 5.8e7
    assert solved["calls"]["n"] == 1, solved["calls"]["n"]
