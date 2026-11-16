"""Example (`EX-124`): ``femem run`` on ``capacitor_sheets`` files -- Route A through the
file equals the circuit reduction, the summary equals its Study, YAML equals Python; then a
file written from the tune summary carries the ruling's 0.2 ohm ESR.

`WF-1` step 13 (item 496) let a YAML file carry ``excitation: {kind: capacitor_sheets, c_f,
esr_ohm, sheets}``: ``run_study_file`` then solves the reduced (1 x 1) network, writes the
run's JSON summary (no Touchstone) and keeps the loaded Study.  No example ran such a file.
This one, on the exact-rectangle PEC hole F-small birdcage (degree 1, 64 MHz, ``-n 2``,
kappa off), through ``cli.main(["run", ...])`` in-process with a call-through wrapper on
``fem_em_solver.study_yaml.run_study_file`` (restored in ``finally``):

* (0) refusals and ``femem check``, before any build: (cX-all), (cX-on) the YAML 1.1 ``on``
  trap, (cX-b1) ``outputs.b1_plus`` refused, (check);
* (1) B1: the tune file -- the 4 x 4 the identity reduces; (R0), (a), (D-disk);
* (2) B2: the tracked capacitor file (lossless, transcribed ``C_mode``): (a'), (I=) the
  summary equals the same Study's result bitwise, (I-A) Route A through the file;
* (3) B3: the Python ``CapacitorSheets`` Study: (I-Py) YAML = Python across two real builds;
* (4) B4: a file this example writes, with ``c_f`` read from the tune run's JSON summary on
  disk and ``esr_ohm`` the ruling's 0.2 ohm: custody by ``float.hex``, run through
  ``femem run``; the ESR readings are PRINTED, never asserted;
* (5) the field of the ESR file's loaded coil (one extra solve on the kept problem):
  (S=-E), (P-E), (cE).

Every record, band and helper is imported from the gate module (the ANS-1 rule); record
digits sit only in the example-local ``REC`` dict with their ``log:line``.  New values print
at <= 6 significant digits; the summary dict, S, Z and z_k are never printed whole.

Negative controls (asserted): (c1) the k = 2 capacitance misses (I-A)'s comparator; (cX-*)
the refusals; (cE) the ESR sheets dissipate (a dropped ``esr_ohm`` would read exactly 0);
(c0) the custody diffs.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:26 -n 2 -t 420

Not claimed: a solved-ESR-through-YAML gate; a tuned f0, Q, R_in, match, efficiency or
absolute S11 / Z_in; copper; degree 2; other frequencies; F-human; a reduced Touchstone; an
external check.
"""

from __future__ import annotations

import contextlib
import copy
import fnmatch
import gc
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
import yaml
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import fem_em_solver.study_yaml as study_yaml  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.lumped import run_lumped_sheet_port_case  # noqa: E402
from fem_em_solver.ports.sparameters import _power_waves  # noqa: E402
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.faraday import (  # noqa: E402
    magnetic_flux_density_from_e,
    project_to_cg1,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

import tests.package.test_wf1_yaml_capacitors as g  # noqa: E402
import tests.validation.test_birdcage_sar_mass_averaged as sar_mod  # noqa: E402
import tests.validation.test_port24_series_rc_sheet as p24  # noqa: E402
import tests.validation.test_port_birdcage_four_port as bc4  # noqa: E402
import tests.validation.test_port_lumped_rlc_termination as rlc  # noqa: E402
from tests.package.test_wf16_study_tune_capacitors import _modal_z, s_to_z  # noqa: E402
from tests.validation.test_birdcage_power_identity import (  # noqa: E402
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_26_femem_run_capacitor_sheets"
RUN_T = OUTPUT_DIR / "mri_26_run_tune"
RUN_C = OUTPUT_DIR / "mri_26_run_capacitors"
RUN_E = OUTPUT_DIR / "mri_26_run_esr"
BAD_DIR = OUTPUT_DIR / "mri_26_bad"
CHAIN_DIR = OUTPUT_DIR / "mri_26_chain"
CHAIN_FILE = CHAIN_DIR / "f_small_hole_pec_capacitors_esr.yaml"
ON_FILE = CHAIN_DIR / "on_key.yaml"

# Gate-log records (docs/testing/logs/20261016T170623Z_WF-1-s13-W1.log), compared at
# 10^-(d-2) with d the printed significant digits (2026-10-14 10:30 ruling's form).
REC = {
    "cells": 80798,                                        # :966-967
    "c1": 6.626796421e-11, "c2": 9.873468900e-11,          # :969 (1e-8)
    "a": -2.976552e-04,                                    # :972 (1e-5)
    "a_prime": 6.139e-11,                                  # :973 (1e-2, magnitude)
    "s11": complex(-2.787323414e-01, 9.590409206e-01),     # :977 (1e-8, |dS|/|S|)
    "z_im": -37.52636761151026,                            # :977 (1e-14); real part 0.0
    "i_a": 1.391811e-05,                                   # :985 (1e-5)
    "c1_factor": 3.5268e-01,                               # :989 (1e-3)
}
R0_PRIOR = (5.052e-11, 3.878e-11, 5.875e-11)              # :971, printed only
I_PY_ONE_MESH = 2.223e-16                                  # :981, printed only
# mri:24 lossless loaded field, EX-120 log :1784 (printed beside, never asserted)
LOSSLESS_CENTRE = {"norm": 1.03902e-06, "b1p": 4.86453e-07}


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def _hex_eq(a, b):
    return float(a).hex() == float(b).hex()


def _refusal(fn):
    """Return the ``StudyConfigError`` ``fn`` raises (exact type), or fail."""
    try:
        fn()
    except g.StudyConfigError as exc:
        assert type(exc) is g.StudyConfigError, type(exc)
        return exc
    raise AssertionError("expected StudyConfigError, nothing raised")


def _run_file(path, outdir, comm):
    """``cli.main(['run', ...])`` with a call-through wrapper; exactly one call asserted."""
    g._fresh_dir(outdir, comm)
    captured = []
    original = study_yaml.run_study_file

    def _call_through(*args, **kwargs):
        out = original(*args, **kwargs)
        captured.append(out)
        return out

    study_yaml.run_study_file = _call_through
    t0 = time.perf_counter()
    try:
        rc = g.cli.main(["run", str(path), "--output-dir", str(outdir)])
    finally:
        study_yaml.run_study_file = original
    wall = time.perf_counter() - t0
    assert study_yaml.run_study_file is original
    rc_all = g._rc_all(rc, comm)
    assert rc_all == (0, 0), rc_all
    assert len(captured) == 1, f"wrapper calls {len(captured)}"
    comm.Barrier()
    return captured[0], wall


def _check_summary(run, outdir, comm):
    """The gate's summary-equals-same-Study body (plain asserts); returns printed facts."""
    s = run.summary
    res = run.study.result
    exc = run.study.excitation
    assert res is run.result
    assert len(s["s_matrix"]) == 1 and len(s["s_matrix"][0]) == 1
    re_im = s["s_matrix"][0][0]
    assert (float(re_im[0]).hex(), float(re_im[1]).hex()) == g._hex(res.s_matrix[0, 0])
    with open(outdir / "summary.json") as fh:
        disk = json.load(fh)
    for key in ("s_matrix", "kept_ports", "capacitors", "cells"):
        assert json.dumps(disk[key]) == json.dumps(s[key]), key
    caps = s["capacitors"]
    assert s["kept_ports"] == ["P1"]
    assert caps["sheets"] == ["P2", "P3", "P4"]
    assert _hex_eq(caps["c_f"], exc.c_f)
    assert _hex_eq(caps["esr_ohm"], exc.esr_ohm)
    zs = complex(caps["impedance_ohm"][0], caps["impedance_ohm"][1])
    z_exc = exc.impedance_ohm(g.FREQUENCY_64_HZ)
    assert g._hex(zs) == g._hex(z_exc) == g._hex(res.capacitor_impedance_ohm)
    assert "touchstone" not in s and "z_matrix" not in s
    listing = sorted((p.name, p.stat().st_size) for p in outdir.iterdir())
    assert not [n for n, _ in listing if fnmatch.fnmatch(n, "*.s*p")], listing
    return listing


def _field_solve(res, ports, z0):
    """Item 373's (P) call (tests/package/test_wf16_study_capacitor_sheets.py:248-251),
    copied from examples/mri/24_study_tune_then_capacitor_sheets.py:116-124."""
    r, f = run_lumped_sheet_port_case(
        res.problem, ports, res.specs, facet_tags=res.facet_tags,
        driven_port_id="P1", verbose=False, return_fields=True,
        extra_bilinear_terms=res.conductor_terms, degree=res.degree)
    resp = r.responses["P1"]
    a, b = _power_waves(resp.voltage_v, resp.current_a, z0)
    return r, f, complex(a), complex(b)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    wall = {}

    # ---- structural asserts -------------------------------------------------
    assert comm.size == g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}")
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    assert g.REDUCTION_BAND is rlc.REDUCTION_BAND
    assert g.LEG_D0_REPRODUCTION_BAND is bc4.LEG_D0_REPRODUCTION_BAND
    assert g.EXACT_IDENTITY_RTOL is sar_mod.EXACT_IDENTITY_RTOL
    assert g.ESR_OHM is p24.ESR_OHM or g.ESR_OHM == p24.ESR_OHM
    band = g.REDUCTION_BAND
    f64 = g.FREQUENCY_64_HZ
    z0 = float(g.REFERENCE_IMPEDANCE_OHM)
    r_esr = float(g.ESR_OHM)

    sy = g.load_study(g.CAP_FILE)
    assert type(sy.excitation) is g.CapacitorSheets
    assert sy.excitation.on == g.SHEETS
    assert sy.excitation.c_f.hex() == g.ITEM429_C_MODE_F.hex()
    assert sy.excitation.esr_ohm == 0.0
    assert sy.degree == 1 and sy.conductor is None and sy.geometry.as_hole is True
    assert g._load(g.CAP_FILE)[2] == {}
    assert g._spec_diff(sy, g._py_caps()) == []
    omega = float(sy.omega)

    # ---- (0) refusals and femem check, before any build ---------------------------
    ex_all = _refusal(lambda: g.load_study(g.MALFORMED_ALL))
    assert ex_all.field == "excitation.sheets", ex_all.field
    OUTPUT_DIR.mkdir(exist_ok=True)
    g._fresh_dir(BAD_DIR, comm)
    with contextlib.redirect_stdout(io.StringIO()):
        rc_bad_run = g.cli.main(["run", str(g.MALFORMED_ALL), "--output-dir", str(BAD_DIR)])
        rc_bad_check = g.cli.main(["check", str(g.MALFORMED_ALL)])
    cx_run, cx_check = g._rc_all(rc_bad_run, comm), g._rc_all(rc_bad_check, comm)
    cx_exists = comm.allreduce(bool(BAD_DIR.exists()), op=MPI.LOR)
    assert cx_run == (2, 2) and cx_check == (2, 2) and not cx_exists, (cx_run, cx_check, cx_exists)

    g._fresh_dir(CHAIN_DIR, comm)  # before (cX-on) writes into it, never after
    if comm.rank == 0:
        CHAIN_DIR.mkdir(parents=True, exist_ok=True)
        ON_FILE.write_text(g.CAP_FILE.read_text().replace(
            "  sheets: [P2, P3, P4]", "  on: [P2, P3, P4]"))
    comm.Barrier()
    ex_on = _refusal(lambda: g.load_study(ON_FILE))
    assert ex_on.field == "excitation.on", ex_on.field
    with contextlib.redirect_stdout(io.StringIO()):
        on_check = g._rc_all(g.cli.main(["check", str(ON_FILE)]), comm)
    assert on_check == (2, 2), on_check

    m_b1 = copy.deepcopy(g._mapping(g.CAP_FILE))
    m_b1["outputs"] = {"b1_plus": {"points": [[0.0, 0.0, 0.0]]}}
    ex_b1 = _refusal(lambda: g._study_from_mapping(m_b1, "<cap + b1_plus>"))
    assert ex_b1.field == "outputs.b1_plus", ex_b1.field

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc_check = g.cli.main(["check", str(g.CAP_FILE)])
    check_all = g._rc_all(rc_check, comm)
    check_lines = buf.getvalue().splitlines()
    assert check_all == (0, 0), check_all
    if comm.rank == 0:
        assert "family: four_port_birdcage" in check_lines
        assert "excitation: CapacitorSheets" in check_lines
    _p0(comm, f"[(0)] cX-all field {ex_all.field!r} rc {cx_run}/{cx_check} dir exists "
              f"{cx_exists}; cX-on field {ex_on.field!r} check rc {on_check}; cX-b1 field "
              f"{ex_b1.field!r}; check rc {check_all}")

    # ---- (1) B1: femem run on the tune file ----------------------------------------
    run_t, wall["B1"] = _run_file(g.TUNE_FILE, RUN_T, comm)
    s0 = np.asarray(run_t.result.s_matrix, dtype=np.complex128).copy()
    assert s0.shape == (4, 4)
    cells_t = int(run_t.summary["cells"])
    assert cells_t == g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD == REC["cells"]
    means = g._class_means(s0)
    want = g.ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]
    r0 = {c: abs(means[c] - want[c]) / abs(want[c]) for c in g.CLASSES}
    c1 = float(run_t.summary["tune"]["c_mode_farad"])
    gap_a = c1 / g.TH17_STEP3_C_MODE_RECORD_F - 1.0
    res_t0 = run_t.study.result
    t2 = g._tune_summary(run_t.study, {"k": 2}, comm)
    assert run_t.study.result is res_t0
    c2 = float(t2["c_mode_farad"])
    with open(RUN_T / "summary.json") as fh:
        disk_t = json.load(fh)
    c_in = float(disk_t["tune"]["c_mode_farad"])
    d_disk = c_in.hex() == c1.hex()
    del run_t, res_t0
    gc.collect()
    comm.Barrier()
    _p0(comm, f"[(R0)] cells {cells_t}; class-mean rel "
              + "  ".join(f"{c} {r0[c]:.3e}" for c in g.CLASSES)
              + f" (band {g.LEG_D0_REPRODUCTION_BAND:.0e}; gate "
              + " / ".join(f"{x:.3e}" for x in R0_PRIOR) + ", :971)")
    _p0(comm, f"[(a)] C_1 = {c1:.9e} F, C_2 = {c2:.9e} F; gap vs TH-17 record {gap_a:+.6e} "
              f"(gate {REC['a']:+.6e}); band {band:.1e}; (D-disk) summary read back == C_1 "
              f"by hex: {d_disk}; B1 {wall['B1']:.1f} s")
    for c in g.CLASSES:
        assert r0[c] <= g.LEG_D0_REPRODUCTION_BAND, f"(R0) {c} {r0[c]:.3e}"
    assert abs(gap_a) <= band, f"(a) {gap_a:+.6e}"
    assert abs(c1 / REC["c1"] - 1.0) <= 1e-8 and abs(c2 / REC["c2"] - 1.0) <= 1e-8
    assert abs(gap_a / REC["a"] - 1.0) <= 1e-5, "(a) record"
    assert d_disk, "(D-disk)"

    # ---- (2) B2: femem run on the tracked capacitor file ---------------------------
    run_c, wall["B2"] = _run_file(g.CAP_FILE, RUN_C, comm)
    resY = run_c.result
    sY = complex(resY.s_matrix[0, 0])
    cells_c = int(run_c.summary["cells"])
    assert cells_c == REC["cells"]
    a_prime = g.ITEM429_C_MODE_F / c1 - 1.0
    listing_c = _check_summary(run_c, RUN_C, comm)
    zc = run_c.study.excitation.impedance_ohm(f64)
    pA = g.reduce_terminated_ports(s0, z0, {1: zc, 2: zc, 3: zc})
    i_a = abs(sY - pA[0, 0]) / abs(pA[0, 0])
    rec_s = abs(sY - REC["s11"]) / abs(REC["s11"])
    rec_zim = abs(zc.imag / REC["z_im"] - 1.0)
    # (c1): the k = 2 capacitance on (I-A)'s comparator
    zc2, pA2 = g._route_a(s0, c2)
    c1_val = abs(sY - pA2[0, 0]) / abs(pA2[0, 0])
    rec_c1 = abs(c1_val / REC["c1_factor"] - 1.0)
    _p0(comm, f"[(a')] ITEM429_C_MODE_F / C_1 - 1 = {a_prime:+.3e} (gate {REC['a_prime']:.3e}); "
              f"inside EXACT_IDENTITY_RTOL {abs(a_prime) <= g.EXACT_IDENTITY_RTOL} (printed); "
              f"band {band:.1e}")
    _p0(comm, f"[(I=)] summary == same Study's result bitwise, JSON read-back, kept ports "
              f"P1, capacitors block, impedance, no Touchstone / z_matrix: all hold; "
              f"directory " + ", ".join(f"{n} {sz} B" for n, sz in listing_c)
              + "; no *.s*p file")
    _p0(comm, f"[S11 file] {sY.real:+.9e}{sY.imag:+.9e}j (rel to :977 {rec_s:.2e}); "
              f"zc = {zc.real:.1f}{zc.imag:+.6f}j ohm (Im rel to :977 {rec_zim:.2e}); "
              f"B2 {wall['B2']:.1f} s")
    _p0(comm, f"[(I-A)] |S_file - reduction|/|reduction| {i_a:.6e} ({i_a / band:.4f} x band "
              f"{band:.0e}); gate :985 {REC['i_a']:.6e}")
    _p0(comm, f"[(c1)] reduction at C_2: {c1_val:.5e} ({c1_val / band:.1f} x band; gate "
              f":989 {REC['c1_factor']:.4e}, rel {rec_c1:.2e}); zc2 != zc {zc2 != zc}")
    assert abs(a_prime) <= band, f"(a') {a_prime:+.6e}"
    assert abs(abs(a_prime) / REC["a_prime"] - 1.0) <= 1e-2, "(a') record"
    assert rec_s <= 1e-8, f"S11 record {rec_s:.2e}"
    assert zc.real == 0.0 and rec_zim <= 1e-14, f"impedance record {zc!r}"
    assert i_a <= band, f"(I-A) {i_a:.6e}"
    assert abs(i_a / REC["i_a"] - 1.0) <= 1e-5, f"(I-A) record {i_a:.6e}"
    assert zc2 != zc and c1_val > band, f"(c1) blind {c1_val:.4e}"
    assert rec_c1 <= 1e-3, f"(c1) record {rec_c1:.2e}"

    # ---- (3) B3: the Python Study ---------------------------------------------------
    py = g._py_caps()
    assert g._spec_diff(run_c.study, py) == []
    t0 = time.perf_counter()
    py_res = py.solve(comm)
    wall["B3"] = time.perf_counter() - t0
    sP = complex(py_res.s_matrix[0, 0])
    i_py = abs(sY - sP) / abs(sP)
    _p0(comm, f"[(I-Py)] YAML vs Python Study (two real builds) |dS|/|S| {i_py:.3e} (band "
              f"{g.EXACT_IDENTITY_RTOL:.0e}; predicted <= 1e-14; gate one-mesh {I_PY_ONE_MESH:.3e}"
              f", :981); hex equal {g._hex(sY) == g._hex(sP)} (printed); spec diff []; "
              f"B3 {wall['B3']:.1f} s")
    assert i_py <= g.EXACT_IDENTITY_RTOL, f"(I-Py) {i_py:.3e}"
    cellsY = resY.cells
    del py, py_res, run_c, resY
    gc.collect()
    comm.Barrier()

    # ---- (4) B4: tune -> place across two files, the ruling's ESR ---------------------
    m = copy.deepcopy(g._mapping(g.CAP_FILE))
    m["excitation"]["c_f"] = c_in
    m["excitation"]["esr_ohm"] = g.ESR_OHM
    if comm.rank == 0:
        CHAIN_FILE.write_text(
            "# mri:26 (EX-124): CAP_FILE with c_f read from the tune run's JSON summary and "
            "the 2026-10-06 ESR ruling's esr_ohm = 0.2 (written by the example).\n"
            + yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    sE = g.load_study(CHAIN_FILE)
    assert sE.excitation.c_f.hex() == c_in.hex()
    assert sE.excitation.esr_ohm.hex() == float(g.ESR_OHM).hex()
    assert sorted(g._spec_diff(sE, sy)) == ["excitation.c_f", "excitation.esr_ohm"]
    m_l = copy.deepcopy(m)
    m_l["excitation"]["esr_ohm"] = 0.0
    s_lossless = g._study_from_mapping(m_l, "<chain lossless>")
    assert g._spec_diff(s_lossless, sy) == ["excitation.c_f"]
    with contextlib.redirect_stdout(io.StringIO()):
        chain_check = g._rc_all(g.cli.main(["check", str(CHAIN_FILE)]), comm)
    assert chain_check == (0, 0), chain_check
    del sE, s_lossless

    run_e, wall["B4"] = _run_file(CHAIN_FILE, RUN_E, comm)
    resE = run_e.result
    sE11 = complex(resE.s_matrix[0, 0])
    listing_e = _check_summary(run_e, RUN_E, comm)
    assert _hex_eq(run_e.summary["capacitors"]["esr_ohm"], g.ESR_OHM)
    zcE = complex(resE.capacitor_impedance_ohm)
    z_ref = g.CapacitorSheets(c_in, g.ESR_OHM, on=g.SHEETS).impedance_ohm(f64)
    assert g._hex(zcE) == g._hex(z_ref)
    assert zcE.real.hex() == float(g.ESR_OHM).hex()
    assert resE.cells == cellsY
    zc0 = g.CapacitorSheets(c_in, 0.0, on=g.SHEETS).impedance_ohm(f64)

    # printed only: the ESR readings
    zs = s_to_z(s0, z0)
    eye = np.eye(4)
    zk0 = np.asarray(_modal_z(zs + zc0 * eye)[0], dtype=np.complex128)
    zkE = np.asarray(_modal_z(zs + zcE * eye)[0], dtype=np.complex128)
    d_re = float(np.max(np.abs(np.real(zkE - zk0) - r_esr)) / r_esr)
    with np.errstate(divide="ignore", invalid="ignore"):
        d_im = float(np.max(np.abs(np.imag(zkE - zk0)) / np.abs(np.imag(zk0))))
    abs_y, abs_e = abs(sY), abs(sE11)
    pAE = g.reduce_terminated_ports(s0, z0, {1: zcE, 2: zcE, 3: zcE})
    i_a_e = abs(sE11 - pAE[0, 0]) / abs(pAE[0, 0])
    _p0(comm, f"[(cE-custody)] written c_f == C_1 by hex, esr_ohm == ruling by hex, spec diff "
              f"vs tracked ['excitation.c_f', 'excitation.esr_ohm'], lossless twin diff "
              f"['excitation.c_f'], check rc {chain_check}")
    _p0(comm, f"[B4 run] (I=) holds on the ESR file; summary esr_ohm == ruling; zc.real == "
              f"ESR by hex; cells == tracked file's; B4 {wall['B4']:.1f} s; directory "
              + ", ".join(f"{n} {sz} B" for n, sz in listing_e))
    _p0(comm, f"[ESR printed] Re lambda shift max_k |Re(zkE - zk0) - r|/r {d_re:.3e} "
              f"(predicted <= 1e-12); max_k |Im(zkE - zk0)|/|Im zk0| {d_im:.3e} (predicted <= "
              f"1e-12)")
    _p0(comm, f"[ESR printed] |S11| lossless {abs_y:.6g}, with ESR {abs_e:.6g} (predicted "
              f"[0.985, 0.997], < lossless: {abs_e < abs_y}); absorbed fraction 1-|S11|^2 "
              f"{1.0 - abs_y ** 2:.4e} -> {1.0 - abs_e ** 2:.4e}; S11 with ESR "
              f"{sE11.real:+.6g}{sE11.imag:+.6g}j")
    _p0(comm, f"[ESR printed] Route A at the lossy termination |S_file - reduction|/|reduction| "
              f"{i_a_e:.4e} ({i_a_e / band:.4f} x band {band:.0e}; predicted [3e-6, 1e-4])")

    # ---- (5) the field of the ESR file's loaded coil ---------------------------------
    assert not resE.conductor_terms
    # setup figure (no-op unless FEM_EM_SETUP_FIGURES=1)
    local_tags = np.unique(resE.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    air_tags = tuple(t for t in all_tags if t != PHANTOM_CELL_TAG)
    names = {PHANTOM_CELL_TAG: "gel phantom"}
    for t in air_tags:
        names[t] = "air (outer box PEC)"
    FIGURE_DIR.mkdir(exist_ok=True)
    write_setup_figure(
        resE.mesh, resE.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=air_tags, translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:26 -- femem run on capacitor_sheets files: C_mode on P2-P4 of the PEC "
               "hole F-small, 0.2 ohm ESR twin (WF-1 step 13)"),
        comm=comm)
    gc.collect()
    comm.Barrier()
    t0 = time.perf_counter()
    rE, fE, aE, bE = _field_solve(resE, resE.all_port_definitions, z0)
    wall["field"] = time.perf_counter() - t0
    s_eq = abs(bE / aE - sE11) / abs(sE11)

    msh = resE.mesh
    e = fE.e_complex
    dx = ufl.Measure("dx", domain=msh, subdomain_data=resE.cell_tags)
    p_vol = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(
        0.5 * fE.sigma_field * ufl.inner(e, e) * dx)), op=MPI.SUM)))
    acct = {"mesh": msh, "facet_tags": resE.facet_tags}
    sheets = {s.port_id: s.sheet(driven=(s.port_id == "P1")) for s in resE.specs}
    p_sheet = {pid: _sheet_field_dissipation_w(acct, sh, e, omega)
               for pid, sh in sheets.items()}
    p_src = float(_source_power_w(acct, sheets["P1"], e, omega))
    p_all = sum(p_sheet.values())
    p_rc = sum(v for pid, v in p_sheet.items() if pid != "P1")
    p_resid = abs(p_src - p_all - p_vol) / abs(p_src)
    share = p_rc / (p_rc + p_vol) if (p_rc + p_vol) != 0.0 else float("nan")

    t0 = time.perf_counter()
    bE_f = project_to_cg1(magnetic_flux_density_from_e(fE.e_complex, omega, name="B_esr_dg0"),
                          name="B_esr")
    bE_f.x.array[:] = bE_f.x.array / aE
    bE_f.x.scatter_forward()
    wall["projection"] = time.perf_counter() - t0
    vals, valid = evaluate_vector_field_parallel(bE_f, np.zeros((1, 3)))
    vals = np.asarray(vals).reshape(-1, 3)
    assert bool(np.all(np.asarray(valid, dtype=bool))), "origin not found"
    v = vals[0]
    centre = {"b1p": abs(v[0] + 1j * v[1]) / 2.0, "b1m": abs(v[0] - 1j * v[1]) / 2.0,
              "norm": float(np.linalg.norm(v))}

    space = bE_f.function_space
    space_s = fem.functionspace(space.mesh, ("Lagrange", 1))

    def _scalar(name, values):
        f = fem.Function(space_s, name=name)
        assert values.shape[0] == f.x.array.shape[0], (values.shape, f.x.array.shape)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    nb = bE_f.x.array.reshape(-1, 3)
    b_real = fem.Function(space, name="B_esr_real")
    b_real.x.array[:] = np.real(bE_f.x.array)
    b_real.x.scatter_forward()
    fields = {
        "B_esr_real": b_real,
        "B_esr_abs": _scalar("B_esr_abs", np.sqrt(np.sum(np.abs(nb) ** 2, axis=1))),
        "B1_plus_esr": _scalar("B1_plus_esr", np.abs(nb[:, 0] + 1j * nb[:, 1]) / 2.0),
        "B1_minus_esr": _scalar("B1_minus_esr", np.abs(nb[:, 0] - 1j * nb[:, 1]) / 2.0),
    }
    for key, fn in fields.items():
        assert fn.name == key, (key, fn.name)
    t0 = time.perf_counter()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, resE.cell_tags, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    wall["xdmf"] = time.perf_counter() - t0
    peak = report_peak_rss(comm, label="EX-124 after outputs")

    # ---- prints (no collective inside) -------------------------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-124 -- femem run on capacitor_sheets files (mri:26)")
        print("=" * 72)
        print(f"[study] degree 1, -n {comm.size}, f = {f64:.6e} Hz, z0 = {z0:g} ohm; tracked "
              f"file lossless at the transcribed C_mode; written file r = {r_esr:g} ohm at "
              f"the in-run C_mode")
        print(f"[(S=-E)] |bE/aE - S11_file|/|S11_file| {s_eq:.3e} (band "
              f"{g.LEG_D0_REPRODUCTION_BAND:.0e}; predicted <= 1e-12; EX-120 lossless 1.842e-13 "
              f":1781), field solve {wall['field']:.1f} s")
        print(f"[(P-E)] residual |P_src - sum P_sheet - P_vol|/P_src {p_resid:.3e} (band "
              f"{DISCRETE_IDENTITY_RTOL:.0e}; predicted <= 1e-12); P_vol > 0 {p_vol > 0.0}")
        print(f"[(cE)] sum P_sheet(P2-P4) > 0: {p_rc > 0.0} (p_rc/P_src {p_rc / p_src:.4e}; "
              f"lossless arm reads exactly 0, EX-120 :1782)")
        print(f"[ESR printed] P_src {p_src:.6e} W; P_sheet(P1)/P_src {p_sheet['P1'] / p_src:.4e}; "
              f"p_rc/P_src {p_rc / p_src:.4e}; P_vol/P_src {p_vol / p_src:.4e}")
        print(f"[ESR printed] ESR share p_rc/(p_rc + P_vol) {share:.4f} (predicted [0.6, 0.95])")
        print(f"[B esr] centre |B| {centre['norm']:.5e}, |B1+| {centre['b1p']:.5e}, "
              f"|B1-| {centre['b1m']:.5e} T per unit incident wave; lossless (mri:24 :1784) "
              f"|B| {LOSSLESS_CENTRE['norm']:.5e}, |B1+| {LOSSLESS_CENTRE['b1p']:.5e}")
        print(f"[B esr] |B_E(0)|/|B_lossless(0)| {centre['norm'] / LOSSLESS_CENTRE['norm']:.4f} "
              f"(predicted [0.98, 1.01])")
        print(f"[wall] B1 {wall['B1']:.1f}, B2 {wall['B2']:.1f}, B3 {wall['B3']:.1f}, B4 "
              f"{wall['B4']:.1f}, field {wall['field']:.1f}, projection "
              f"{wall['projection']:.1f}, xdmf {wall['xdmf']:.1f} s; summed peak RSS "
              f"{peak / 2**30:.2f} GiB")
        print(f"[paraview] wrote {xdmf_path}")

    # ---- asserts for the field arm -------------------------------------------------------
    assert s_eq <= g.LEG_D0_REPRODUCTION_BAND, f"(S=-E) {s_eq:.3e}"
    assert p_vol > 0.0
    assert p_resid <= DISCRETE_IDENTITY_RTOL, f"(P-E) {p_resid:.3e}"
    assert p_rc > 0.0, "(cE) the ESR sheets dissipate nothing: esr_ohm did not reach the field"

    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
