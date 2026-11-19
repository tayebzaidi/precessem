"""`WF-1` step 12: a top-level ``tune:`` on the copper F-human hole file drives `WF-16`'s ring read.

``run_study_file`` on ``f_human_hole_copper_tune.yaml`` (``f_human_birdcage``
hole, ``lumped_port_drive``, Leontovich copper, degree 1, 64 MHz, κ off,
``-n 8``) solves the 32-port ring-gap matrix once and writes ``Study.tune(k=1)``
into the summary; its ``C_mode`` is compared with `WF-12` step 3's H12c reading
within ``REDUCTION_BAND`` (brief: docs/planning/chunks/WF-1.md, §9 item 484).

Anchors (asserted, every band imported, none set here):
  (K)   the parse equals item 472's Python Study field for field; no mesh.
  (X)   every refusal names its field; no mesh.
  (R0)  cells == ``HOLE_CELL_RECORD``; the run's 32x32 ties the H12c npz at
        ``REUSE_REPRODUCTION_RTOL``; a ``.s32p`` is written.
  (I=)  the summary equals the same Study's ``tune(k=1)`` at
        ``EXACT_IDENTITY_RTOL`` (one solve; ``float.hex`` printed).
  (a)   ``C_mode`` against the H12c record within ``REDUCTION_BAND``.
Negative control (rule (e)), asserted:
  (c1)  k = 2 through the YAML route's summary function misses the record.

Opt-in: ``FEM_EM_WF1_S12=1`` (unset => the solve nodes skip; any other value
raises); width 8 asserted when set.  Scope: degree 1, 64 MHz, κ off, lossless
read, one copper hole fixture; not a tuned f0, Q, R_in, C_in or resonance, not
ESR, not degree 2, not AED.
"""

from __future__ import annotations

import copy
import io
import json
import os
import time
from contextlib import redirect_stdout

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import StudyConfigError, cli, load_study, run_study_file, study_yaml
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    SphereInBoxGeometry,
)
from fem_em_solver.study_yaml import _parse_mapping, _study_from_mapping, _tune_summary
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf16_study_tune_f_human_copper import (
    RECORD_RANK_WIDTH,
    WF12_STEP3_C_MODE_CU_F,
    WF12_STEP3_K2_SEPARATION,
    _study,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_wf10_f_human_copper_esr_root import NPZ_PATH
from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

S12_ENV = "FEM_EM_WF1_S12"
FILE = WF1 / "f_human_hole_copper_tune.yaml"
QUAD_FILE = WF1 / "f_human_tuned_copper.yaml"
#: Item 472's measured API reading on this fixture and width, printed beside only:
#: C_mode and gap (20261015T213534Z_WF-16-s12-W1.log:7375), k = 2 gap (:7380).
ITEM472_C_MODE_F = 1.181185799070e-10
ITEM472_GAP = -3.544942e-13
ITEM472_K2_GAP = -0.3327259


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S12_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S12_ENV}={raw!r}: only unset/empty or '1' are accepted")


def test_wf1_s12_parse_and_refusals(monkeypatch):
    """(K) and (X): no mesh, no solve."""
    _knob()  # a bad knob value raises here too
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)

    # (K)
    d_py = _spec_diff(load_study(FILE), _study())
    out = study_yaml._load(FILE)[2]
    _p0(f"\n[WF-1 s12] (K) diff vs item 472's _study() {d_py}; parsed outputs {out!r}")
    assert d_py == [], d_py
    assert out == {"tune": {"k": 1}}
    assert type(out["tune"]["k"]) is int
    m_empty = copy.deepcopy(_mapping(FILE))
    m_empty["tune"] = {}
    assert _parse_mapping(m_empty, "<tune {}>")[2] == {"tune": {}}

    # (X) mapping variants of FILE, each naming its field
    base = _mapping(FILE)

    def variant(fn):
        m = copy.deepcopy(base)
        fn(m)
        return m

    def _add_thermal(m):  # item 390a's values (test_wf1_yaml_f_human.py:215–216)
        m["thermal"] = {"k_w_per_m_k": 0.6, "h_w_per_m2_k": 10.0, "t_inf_c": 20.0}

    cases = [
        ("tune k: 3", variant(lambda m: m.__setitem__("tune", {"k": 3})), "tune.k"),
        (
            "excitation.c_in_f added",
            variant(lambda m: m["excitation"].__setitem__("c_in_f", 1.18e-10)),
            "excitation.kind",
        ),
        ("thermal added", variant(_add_thermal), "thermal"),
        (
            "outputs.sar",
            variant(
                lambda m: m.__setitem__(
                    "outputs", {"sar": {"mass_g": 10.0, "centres": [[0.0147, 0.0029, 0.0]]}}
                )
            ),
            "outputs.sar",
        ),
        (
            "outputs.b1_plus",
            variant(lambda m: m.__setitem__("outputs", {"b1_plus": {"points": [[0.0, 0.0, 0.0]]}})),
            "outputs.b1_plus",
        ),
        (
            "as_hole false",
            variant(lambda m: m["geometry"].__setitem__("as_hole", False)),
            "geometry.as_hole",
        ),
        ("conductor removed", variant(lambda m: m.pop("conductor")), "conductor"),
    ]
    m_quad = copy.deepcopy(_mapping(QUAD_FILE))
    m_quad["tune"] = {"k": 1}
    cases.append(("tune on f_human_tuned_copper.yaml", m_quad, "tune"))
    fields, expected = [], []
    for label, m, want in cases:
        with pytest.raises(StudyConfigError) as ei:
            _study_from_mapping(m, f"<{label}>")
        _p0(f"  (X) {label}: field {ei.value.field!r} (expected {want!r}); {ei.value}")
        assert type(ei.value) is StudyConfigError
        fields.append(ei.value.field)
        expected.append(want)
    assert fields == expected, fields

    # cli check on the good file: rc 0 on every rank, family printed on rank 0
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = cli.main(["check", str(FILE)])
    lo_hi = _rc_all(rc, comm)
    text = buf.getvalue()
    _p0(f"  check rc (min, max) {lo_hi}; rank-0 stdout {text.strip()!r}")
    assert lo_hi == (0, 0)
    if comm.rank == 0:
        assert "family: f_human_birdcage" in text
    assert calls["n"] == 0


@pytest.fixture(scope="module")
def runs():
    if not _knob():
        pytest.skip(
            f"{S12_ENV} unset: the copper F-human YAML tune is a heavy -n 8 window, "
            "environment-selected"
        )
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, (comm.size, RECORD_RANK_WIDTH)
    npz = None
    if comm.rank == 0:
        assert os.path.exists(NPZ_PATH), f"{NPZ_PATH} absent: stop (BLOCKED, rule (d))"
        with np.load(NPZ_PATH, allow_pickle=False) as fh:
            npz = {
                "s32": np.array(fh["s32"], dtype=np.complex128),
                "order": [str(x) for x in fh["order"]],
                "cells": int(fh["cells"]),
            }
    npz = comm.bcast(npz, root=0)
    d = LOGS / "f_human_hole_copper_tune"
    _fresh_dir(d, comm)
    t0 = time.perf_counter()
    run = run_study_file(FILE, output_dir=d, comm=comm)
    t_run = time.perf_counter() - t0
    comm.Barrier()
    with open(d / "summary.json") as fh:
        summary = json.load(fh)
    res0 = run.study.result
    api1 = run.study.tune(k=1, comm=comm)
    assert run.study.result is res0  # one solve
    ts_exists = comm.bcast(
        bool((d / summary["touchstone"]).exists()) if comm.rank == 0 else None, root=0
    )
    _p0(
        f"\n[WF-1 s12] cells {res0.cells}, mesh {res0.mesh_time:.1f} s, sweep "
        f"{res0.sweep_time:.1f} s; run_study_file {t_run:.1f} s; -n {comm.size}"
    )
    report_peak_rss(comm, label="[WF-1 s12] summed peak ru_maxrss")
    _p0(f"[WF-1 s12] summary tune {summary['tune']}")
    return {
        "run": run,
        "summary": summary,
        "res0": res0,
        "api1": api1,
        "npz": npz,
        "ts_exists": ts_exists,
    }


@complex_only
def test_wf1_s12_r0_fixture_tie(runs):
    """(R0): the YAML Study is the H12c fixture; its 32x32 ties the stored matrix."""
    run, summary, npz = runs["run"], runs["summary"], runs["npz"]
    s = np.asarray(run.result.s_matrix, dtype=np.complex128)
    names = [pd.port_id for pd in run.result.port_definitions]
    idx = {n: i for i, n in enumerate(names)}
    p = np.array([idx[n] for n in npz["order"]], dtype=int)
    s_p = s[np.ix_(p, p)]
    rel = _residual(s_p, npz["s32"])
    col = np.linalg.norm(s_p - npz["s32"], axis=0) / np.linalg.norm(npz["s32"], axis=0)
    _p0(
        f"\n[WF-1 s12 R0] cells {summary['cells']} vs record {HOLE_CELL_RECORD} vs npz "
        f"{npz['cells']}; ||S_run - S_npz||/||S_npz|| = {rel:.3e} (band "
        f"REUSE_REPRODUCTION_RTOL {REUSE_REPRODUCTION_RTOL:.1e}; predicted <= 1e-10); worst "
        f"column {col.max():.3e}; touchstone {summary['touchstone']} exists {runs['ts_exists']}"
    )
    assert summary["cells"] == HOLE_CELL_RECORD
    assert sorted(names) == sorted(npz["order"])
    assert rel <= REUSE_REPRODUCTION_RTOL
    assert len(summary["s_matrix"]) == 32
    assert summary["touchstone"] == FILE.stem + ".s32p"
    assert runs["ts_exists"]


@complex_only
def test_wf1_s12_summary_equals_same_study_tune(runs):
    """(I=): the summary is the same Study's ``tune(k=1)``; one solve."""
    run, summary, api1 = runs["run"], runs["summary"], runs["api1"]
    mem = run.summary["tune"]
    rc = abs(mem["c_mode_farad"] / api1.c_mode_farad - 1.0)
    ro = abs(mem["off_diagonal_residual"] / api1.off_diagonal_residual - 1.0)
    _p0(
        f"\n[WF-1 s12 I=] summary C {mem['c_mode_farad']!r} vs api {api1.c_mode_farad!r}: rel "
        f"{rc:.3e}, off rel {ro:.3e} (band EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.1e}); hex "
        f"{mem['c_mode_farad'].hex()} vs {float(api1.c_mode_farad).hex()} equal "
        f"{mem['c_mode_farad'].hex() == float(api1.c_mode_farad).hex()} (printed)"
    )
    assert rc <= EXACT_IDENTITY_RTOL
    assert ro <= EXACT_IDENTITY_RTOL
    assert len(mem["z_k"]) == 9 and len(api1.z_k) == 9
    for (re, im), z in zip(mem["z_k"], api1.z_k):
        z = complex(z)
        assert abs(re - z.real) <= EXACT_IDENTITY_RTOL * abs(z)
        assert abs(im - z.imag) <= EXACT_IDENTITY_RTOL * abs(z)
    assert summary["tune"] == mem  # JSON read-back
    assert mem["k"] == 1
    assert mem["frequency_hz"] == 64e6


@complex_only
def test_wf1_s12_a_c_mode_equals_wf12_h12c(runs):
    """(a): the YAML C_mode within REDUCTION_BAND of `WF-12` step 3's H12c reading."""
    c = runs["summary"]["tune"]["c_mode_farad"]
    gap = c / WF12_STEP3_C_MODE_CU_F - 1.0
    _p0(
        f"\n[WF-1 s12 a] C_mode = {c:.12e} F; record {WF12_STEP3_C_MODE_CU_F:.12e} F; gap "
        f"{gap:+.6e} (band REDUCTION_BAND {REDUCTION_BAND:.1e}; predicted <= 1e-9); vs item "
        f"472's C {c / ITEM472_C_MODE_F - 1:+.3e} (item 472 gap {ITEM472_GAP:+.6e}, printed)"
    )
    assert abs(gap) <= REDUCTION_BAND


@complex_only
def test_wf1_s12_c1_mode_index_control(runs):
    """(c1): ``k = 2`` through ``_tune_summary`` misses the record (predicted -33.27 %)."""
    comm = MPI.COMM_WORLD
    t2 = _tune_summary(runs["run"].study, {"k": 2}, comm)
    assert runs["run"].study.result is runs["res0"]
    assert t2["k"] == 2
    gap = t2["c_mode_farad"] / WF12_STEP3_C_MODE_CU_F - 1.0
    _p0(
        f"\n[WF-1 s12 c1] k = 2 C = {t2['c_mode_farad']:.12e} F, gap {gap:+.6e} "
        f"({abs(gap) / REDUCTION_BAND:.1f}x band; item 472 {ITEM472_K2_GAP:+.7f}; step 3 "
        f"separation {WF12_STEP3_K2_SEPARATION})"
    )
    assert abs(gap) > REDUCTION_BAND
