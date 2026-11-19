"""`WF-1` step 11: a top-level ``tune:`` section drives `WF-16`'s ``Study.tune``.

``run_study_file`` on the F-small exact-rectangle PEC hole Study (degree 1,
64 MHz, ``-n 2``, κ off) calls ``Study.tune(k)`` on the solved Study and writes
the result into the summary; its ``C_mode`` is compared with `TH-17`'s step-3
record within ``REDUCTION_BAND`` (brief: docs/planning/chunks/WF-1.md).

Anchors: (K) the parse is bitwise and needs no mesh; (X) every refusal names
its field and builds no mesh; (R0) the YAML Study is the registered hole
fixture; (I=) the summary equals the same Study's ``tune()`` by ``float.hex``;
(I-Py) it equals an independently solved Python Study's within
``LEG_D0_REPRODUCTION_BAND`` (rule (j): bitwise stops at the solve); (a) the
record comparison.  Negative control (asserted): (c1) ``k = 2`` through the
YAML route's summary function misses the record by more than the band.  Every
band and record is imported.
"""

from __future__ import annotations

import copy
import json
import time

import pytest
from mpi4py import MPI

from fem_em_solver import StudyConfigError, cli, load_study, run_study_file, study_yaml
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    SphereInBoxGeometry,
)
from fem_em_solver.study_yaml import _parse_mapping, _study_from_mapping, _tune_summary
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_copper import PEC_CONTROL_FILE, _py_copper
from tests.package.test_wf16_study_copper import _rel_means
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

KEY = "64"
TUNE_FILE = WF1 / "f_small_hole_pec_tune.yaml"
MALFORMED_TUNE = WF1 / "malformed_tune_k_zero.yaml"
# Item 429's measured (a) on this fixture, degree and width (printed beside, never asserted).
ITEM429_A_PRINT = -2.976552e-4


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _hex(x):
    return float(x).hex()


def test_wf1_s11_tune_parse_and_refusals(monkeypatch):
    """(K), (X), (c0): no mesh, no solve."""
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)

    # (K)
    d_py = _spec_diff(load_study(TUNE_FILE), _py_copper("pec"))
    d_ctl = _spec_diff(load_study(TUNE_FILE), load_study(PEC_CONTROL_FILE))
    out_tune = study_yaml._load(TUNE_FILE)[2]
    out_ctl = study_yaml._load(PEC_CONTROL_FILE)[2]
    _p0(
        f"\n[WF-1 s11] (K) diff vs _py_copper('pec') {d_py}; vs PEC control {d_ctl}; "
        f"parsed outputs {out_tune!r}; control outputs {out_ctl!r}"
    )
    assert d_py == [], d_py
    assert d_ctl == [], d_ctl
    assert out_tune == {"tune": {"k": 1}}
    assert type(out_tune["tune"]["k"]) is int
    assert "tune" not in out_ctl  # (c0)
    m_empty = copy.deepcopy(_mapping(TUNE_FILE))
    m_empty["tune"] = {}
    assert _parse_mapping(m_empty, "<tune {}>")[2] == {"tune": {}}
    m_k2 = copy.deepcopy(_mapping(TUNE_FILE))
    m_k2["tune"] = {"k": 2}
    assert _parse_mapping(m_k2, "<tune k2>")[2] == {"tune": {"k": 2}}

    # (X) the tracked malformed file, through load_study and cli.main
    with pytest.raises(StudyConfigError) as ei:
        load_study(MALFORMED_TUNE)
    _p0(f"  (X) {MALFORMED_TUNE.name}: field {ei.value.field!r}; {ei.value}")
    assert type(ei.value) is StudyConfigError
    assert ei.value.field == "tune.k"
    out = LOGS / "malformed_tune"
    _fresh_dir(out, comm)
    rc_run = cli.main(["run", str(MALFORMED_TUNE), "--output-dir", str(out)])
    rc_check = cli.main(["check", str(MALFORMED_TUNE)])
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
    _p0(f"  (X) run rc (min, max) {run_all}, check rc {check_all}, out exists {exists}")
    assert run_all == (2, 2)
    assert check_all == (2, 2)
    assert not exists

    # (X) mapping variants
    base = _mapping(TUNE_FILE)
    variants = []
    for label, tune_value in (
        ("k: 3", {"k": 3}),
        ("k: 1.0", {"k": 1.0}),
        ("k: true", {"k": True}),
        ("mode: 1", {"mode": 1}),
        ("tune: 1", 1),
    ):
        m = copy.deepcopy(base)
        m["tune"] = tune_value
        variants.append((label, m))
    for name in ("sphere_mat4.yaml", "f_human_tuned_copper.yaml"):
        m = copy.deepcopy(_mapping(WF1 / name))
        m["tune"] = {"k": 1}
        variants.append((f"tune on {name}", m))
    expected = ["tune.k", "tune.k", "tune.k", "tune.mode", "tune", "tune", "tune"]
    fields = []
    for (label, m), want in zip(variants, expected):
        with pytest.raises(StudyConfigError) as ei:
            _study_from_mapping(m, f"<{label}>")
        _p0(f"  (X) {label}: field {ei.value.field!r} (expected {want!r}); {ei.value}")
        assert type(ei.value) is StudyConfigError
        fields.append(ei.value.field)
    assert fields == expected, fields

    m_top = copy.deepcopy(base)
    m_top["tunes"] = {"k": 1}
    with pytest.raises(StudyConfigError) as ei:
        _study_from_mapping(m_top, "<tunes>")
    _p0(f"  (X) unknown top-level key message: {ei.value}")
    assert ei.value.field == "tunes"
    assert calls["n"] == 0


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    d = LOGS / "f_small_hole_pec_tune"
    _fresh_dir(d, comm)
    t0 = time.perf_counter()
    run = run_study_file(TUNE_FILE, output_dir=d, comm=comm)
    t_run = time.perf_counter() - t0
    with open(d / "summary.json") as fh:
        summary = json.load(fh)
    res0 = run.study.result
    api1 = run.study.tune(k=1, comm=comm)
    same_after_api = run.study.result is res0
    t0 = time.perf_counter()
    py = _py_copper("pec")
    py_res = py.solve(comm)
    t_py = time.perf_counter() - t0
    pyt = py.tune(k=1, comm=comm)
    _p0(
        f"\n[WF-1 s11] YAML build: cells {run.result.cells}, mesh {run.result.mesh_time:.1f} s, "
        f"sweep {run.result.sweep_time:.1f} s; run_study_file {t_run:.1f} s; "
        f"summary elapsed_s {summary.get('elapsed_s')}"
    )
    _p0(
        f"[WF-1 s11] Python build: cells {py_res.cells}, mesh {py_res.mesh_time:.1f} s, "
        f"sweep {py_res.sweep_time:.1f} s; total {t_py:.1f} s"
    )
    _p0(f"[WF-1 s11] summary tune {summary['tune']}")
    return {
        "run": run,
        "summary": summary,
        "res0": res0,
        "api1": api1,
        "same_after_api": same_after_api,
        "py": py,
        "pyt": pyt,
    }


@complex_only
def test_wf1_s11_r0_fixture_tie(runs):
    """(R0): the YAML Study is the registered exact-rectangle PEC hole fixture."""
    run, summary = runs["run"], runs["summary"]
    assert summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    rel = _rel_means(_class_means(run.result.s_matrix), ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["pec"])
    _p0(
        f"\n[WF-1 s11 R0] cells {summary['cells']}; class-mean rel vs 'pec' record: "
        + ", ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])


@complex_only
def test_wf1_s11_yaml_equals_same_study_api(runs):
    """(I=): the summary is the same Study's ``tune()``, bitwise; no second solve."""
    run, summary, api1 = runs["run"], runs["summary"], runs["api1"]
    assert runs["same_after_api"]
    mem = run.summary["tune"]
    _p0(f"\n[WF-1 s11 I=] summary C {mem['c_mode_farad'].hex()} vs api {api1.c_mode_farad.hex()}")
    assert mem["c_mode_farad"].hex() == api1.c_mode_farad.hex()
    assert mem["off_diagonal_residual"].hex() == api1.off_diagonal_residual.hex()
    assert len(mem["z_k"]) == len(api1.z_k)
    for (re, im), z in zip(mem["z_k"], api1.z_k):
        assert re.hex() == _hex(complex(z).real)
        assert im.hex() == _hex(complex(z).imag)
    js = summary["tune"]
    assert js["c_mode_farad"].hex() == mem["c_mode_farad"].hex()
    assert js["off_diagonal_residual"].hex() == mem["off_diagonal_residual"].hex()
    assert js["frequency_hz"].hex() == mem["frequency_hz"].hex()
    assert js["k"] == mem["k"]
    for (a, b), (c, e) in zip(js["z_k"], mem["z_k"]):
        assert a.hex() == c.hex() and b.hex() == e.hex()
    assert js["k"] == 1
    assert float(js["frequency_hz"]).hex() == FREQUENCY_64_HZ.hex()


@complex_only
def test_wf1_s11_yaml_equals_python_route(runs):
    """(I-Py): YAML vs an independent Python Study, across builds, within the band."""
    c_yaml = runs["summary"]["tune"]["c_mode_farad"]
    c_py = runs["pyt"].c_mode_farad
    rel = abs(c_yaml / c_py - 1.0)
    _p0(
        f"\n[WF-1 s11 I-Py] YAML C {c_yaml:.12e} vs Python C {c_py:.12e}: rel {rel:.3e} "
        f"(band {LEG_D0_REPRODUCTION_BAND:.0e}; predicted <= 1e-14); hex equal "
        f"{c_yaml.hex() == c_py.hex()} (printed, never asserted)"
    )
    assert rel <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf1_s11_c_mode_equals_th17_record(runs):
    """(a): the YAML C_mode within REDUCTION_BAND of `TH-17`'s record."""
    c = runs["summary"]["tune"]["c_mode_farad"]
    gap = c / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(
        f"\n[WF-1 s11 a] C_mode = {c:.9e} F; record {TH17_STEP3_C_MODE_RECORD_F:.9e} F; gap "
        f"{gap:+.6e} (predicted -2.977e-4 in [-3.0e-4, -2.95e-4]; item 429 measured "
        f"{ITEM429_A_PRINT:+.6e}, printed); band {REDUCTION_BAND:.1e}"
    )
    assert abs(gap) <= REDUCTION_BAND


@complex_only
def test_wf1_s11_c1_mode_index_control(runs):
    """(c1): ``k = 2`` through the YAML route's summary function misses the record."""
    comm = MPI.COMM_WORLD
    m = copy.deepcopy(_mapping(TUNE_FILE))
    m["tune"] = {"k": 2}
    spec = _parse_mapping(m, "<tune k2>")[2]["tune"]
    assert spec == {"k": 2}
    t2 = _tune_summary(runs["run"].study, spec, comm)
    assert runs["run"].study.result is runs["res0"]
    assert t2["k"] == 2
    gap = t2["c_mode_farad"] / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(
        f"\n[WF-1 s11 c1] k = 2 C = {t2['c_mode_farad']:.9e} F, gap {gap:+.6e} "
        f"({abs(gap) / REDUCTION_BAND:.1f}x band; item 429 measured +48.95 %)"
    )
    assert abs(gap) > REDUCTION_BAND
