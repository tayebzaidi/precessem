"""`WF-1` step 13: a ``capacitor_sheets`` excitation drives `WF-16`'s ``CapacitorSheets`` from YAML.

``run_study_file`` on the F-small exact-rectangle PEC hole Study (degree 1,
64 MHz, ``-n 2``, κ off, no conductor) with lossless capacitor sheets on P2-P4
at the k = 1 ``C_mode`` reproduces the Route-A circuit = field identity
(`PORT-22`'s construction) against the tune file's own untuned 4x4 within
``REDUCTION_BAND``; its summary equals the same Study's result by
``float.hex`` and the Python ``CapacitorSheets`` Study on the same mesh within
``EXACT_IDENTITY_RTOL``.  The tracked ``c_f`` is a transcribed record (the
10-digit print of the in-run `TH-17` / ``Study.tune(k=1)`` ``C_mode``) whose
custody (a') checks against this run's own tune.  The build cache (one real
mesh, later equal-geometry builds return it) is a test economy licensed by
(K)'s spec equality; its counts are asserted.  Every band and record is
imported; none is set here (brief: docs/planning/chunks/WF-1.md).

Anchors: (K) parse, no mesh; (X) every refusal names its field and builds no
mesh; (R0) the fixture tie and (a) / (a'); (I=) summary = the same Study's
result, bitwise; (I-Py) YAML Study = Python Study at ``EXACT_IDENTITY_RTOL``;
(I-A) the Route-A reduction at ``REDUCTION_BAND``.  Negative controls: (c1)
the k = 2 capacitance on (I-A)'s comparator misses the band (asserted); (c0)
the dropped-key diffs (asserted, by construction).
"""

from __future__ import annotations

import copy
import json
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    Study,
    StudyConfigError,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import reduce_terminated_ports
from fem_em_solver.study import FHumanBirdcageGeometry, SphereInBoxGeometry
from fem_em_solver.study_yaml import _load, _study_from_mapping, _tune_summary
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_copper import PEC_CONTROL_FILE
from tests.package.test_wf1_yaml_tune import TUNE_FILE
from tests.package.test_wf16_study_capacitor_sheets import C_F, ESR_OHM, _study
from tests.package.test_wf16_study_copper import _geometry
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
from tests.package.test_wf16_study_tune_capacitors import ITEM429_C_MODE_F
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

CAP_FILE = WF1 / "f_small_hole_pec_capacitors.yaml"
MALFORMED_ALL = WF1 / "malformed_capacitor_sheets_all_ports.yaml"
SHEETS = ("P2", "P3", "P4")
# Item 463's API-route readings on this fixture and width (printed beside, never asserted).
ITEM463_IA = 1.391811e-05  # 20261015T201110Z_WF-16-s6b-W1.log:954
ITEM463_S11 = complex(-2.787323415e-01, 9.590409206e-01)  # same line
ITEM463_C2 = 3.5268e-01  # :971
ITEM463_A_PRIME = -6.139e-11  # :950


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _say(capsys, *args):
    # past capsys, so the record lands in the harness log
    with capsys.disabled():
        _p0(*args)


def _hex(z):
    z = complex(z)
    return (z.real.hex(), z.imag.hex())


def _py_caps(c_f=ITEM429_C_MODE_F, esr_ohm=0.0):
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    return Study(
        _geometry(),
        phantom,
        CapacitorSheets(c_f, esr_ohm, on=SHEETS),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
    )


def _with_exc(path, **exc):
    m = copy.deepcopy(_mapping(path))
    m["excitation"] = {**m["excitation"], **exc}
    return m


def test_wf1_s13_parse_and_refusals(monkeypatch, capsys):
    """(K), (X), (c0): no mesh, no solve."""
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)

    # (K)
    s_yaml = load_study(CAP_FILE)
    d_py = _spec_diff(s_yaml, _py_caps())
    d_ctl = _spec_diff(s_yaml, load_study(PEC_CONTROL_FILE))  # (c0): a dropped section
    _say(
        capsys,
        f"\n[WF-1 s13] (K) diff vs Python CapacitorSheets Study {d_py}; vs PEC "
        f"control {d_ctl}; c_f {s_yaml.excitation.c_f.hex()} vs ITEM429_C_MODE_F "
        f"{ITEM429_C_MODE_F.hex()}; degree {s_yaml.degree}; conductor {s_yaml.conductor}",
    )
    assert d_py == [], d_py
    assert type(s_yaml.excitation) is CapacitorSheets
    assert s_yaml.excitation.on == SHEETS
    assert s_yaml.excitation.c_f.hex() == ITEM429_C_MODE_F.hex()
    assert s_yaml.degree == 1 and s_yaml.conductor is None
    assert d_ctl == ["excitation"], d_ctl
    assert _load(CAP_FILE)[2] == {}

    # (K-esr) and (c0): ESR reaches the Study
    s_esr = _study_from_mapping(_with_exc(CAP_FILE, esr_ohm=0.2), "<esr 0.2>")
    d_esr = _spec_diff(s_esr, s_yaml)
    z_yaml = s_esr.excitation.impedance_ohm(FREQUENCY_64_HZ)
    z_api = CapacitorSheets(ITEM429_C_MODE_F, 0.2, on=SHEETS).impedance_ohm(FREQUENCY_64_HZ)
    _say(
        capsys, f"  (K-esr) diff {d_esr}; impedance yaml {complex(z_yaml)!r} api {complex(z_api)!r}"
    )
    assert d_esr == ["excitation.esr_ohm"], d_esr
    assert _hex(z_yaml) == _hex(z_api)

    # (K-cu): item 373's copper ESR Study expressed from YAML, parse only
    m_cu = copy.deepcopy(_mapping(WF1 / "f_small_copper.yaml"))
    m_cu["excitation"] = {
        "kind": "capacitor_sheets",
        "c_f": C_F,
        "esr_ohm": ESR_OHM,
        "sheets": ["P2"],
    }
    d_cu = _spec_diff(
        _study_from_mapping(m_cu, "<cu>"), _study(CapacitorSheets(C_F, ESR_OHM, on=("P2",)))
    )
    _say(
        capsys,
        f"  (K-cu) diff vs item 373's copper Study {d_cu} (C_F {C_F!r}, ESR_OHM {ESR_OHM!r})",
    )
    assert d_cu == [], d_cu

    # omitted esr_ohm: the dataclass default, not filled by the parser
    m_no = copy.deepcopy(_mapping(CAP_FILE))
    del m_no["excitation"]["esr_ohm"]
    s_no = _study_from_mapping(m_no, "<no esr>")
    assert s_no.excitation.esr_ohm == 0.0
    assert _spec_diff(s_no, s_yaml) == []

    # cli check
    capsys.readouterr()
    rc = cli.main(["check", str(CAP_FILE)])
    out = capsys.readouterr().out.splitlines()
    rc_all = _rc_all(rc, comm)
    _say(capsys, f"  (K) femem check rc {rc_all}; rank-0 lines {out}")
    assert rc_all == (0, 0)
    if comm.rank == 0:
        assert "family: four_port_birdcage" in out
        assert "excitation: CapacitorSheets" in out

    # (X) the tracked malformed file, through load_study and cli.main
    with pytest.raises(StudyConfigError) as ei:
        load_study(MALFORMED_ALL)
    _say(capsys, f"  (X) {MALFORMED_ALL.name}: field {ei.value.field!r}; {ei.value}")
    assert type(ei.value) is StudyConfigError
    assert ei.value.field == "excitation.sheets"
    d = LOGS / "malformed_capacitor_sheets"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED_ALL), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED_ALL)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(capsys, f"  (X) run rc {run_all}, check rc {check_all}, out exists {exists}")
    assert run_all == (2, 2)
    assert check_all == (2, 2)
    assert not exists

    # (X) the YAML 1.1 `on` key through a real YAML text
    scratch = LOGS / "s13_on_key"
    _fresh_dir(scratch, comm)
    on_file = scratch / "on_key.yaml"
    if comm.rank == 0:
        scratch.mkdir(parents=True, exist_ok=True)
        on_file.write_text(
            CAP_FILE.read_text().replace("  sheets: [P2, P3, P4]", "  on: [P2, P3, P4]")
        )
    comm.Barrier()
    with pytest.raises(StudyConfigError) as ei:
        load_study(on_file)
    _say(capsys, f"  (X) `on:` YAML text: field {ei.value.field!r}; {ei.value}")
    assert ei.value.field == "excitation.on"

    # (X) mapping variants
    base = _mapping(CAP_FILE)
    variants = []

    def _v(label, m, want):
        variants.append((label, m, want))

    for label, value, want in (
        ("sheets: [P5]", ["P5"], "excitation.sheets[0]"),
        ("sheets: [P2, P2]", ["P2", "P2"], "excitation.sheets"),
        ("sheets: []", [], "excitation.sheets"),
        ("sheets: P2", "P2", "excitation.sheets"),
    ):
        _v(label, _with_exc(CAP_FILE, sheets=value), want)
    m = copy.deepcopy(base)
    del m["excitation"]["sheets"]
    _v("sheets removed", m, "excitation.sheets")
    m = copy.deepcopy(base)
    m["excitation"][True] = m["excitation"].pop("sheets")
    _v("True key", m, "excitation.on")
    _v("c_f: 0.0", _with_exc(CAP_FILE, c_f=0.0), "excitation.c_f")
    _v("c_f: '66e-12'", _with_exc(CAP_FILE, c_f="66e-12"), "excitation.c_f")
    m = copy.deepcopy(base)
    del m["excitation"]["c_f"]
    _v("c_f removed", m, "excitation.c_f")
    _v("esr_ohm: -0.1", _with_exc(CAP_FILE, esr_ohm=-0.1), "excitation.esr_ohm")
    _v("z0_ohm added", _with_exc(CAP_FILE, z0_ohm=50.0), "excitation.z0_ohm")
    m = copy.deepcopy(base)
    m["geometry"] = copy.deepcopy(_mapping(WF1 / "sphere_mat4.yaml")["geometry"])
    _v("sphere geometry", m, "excitation.kind")
    m = copy.deepcopy(_mapping(WF1 / "f_human_hole_copper_tune.yaml"))
    m.pop("tune", None)
    m["excitation"] = copy.deepcopy(base["excitation"])
    _v("f_human geometry", m, "excitation.kind")
    m = copy.deepcopy(base)
    m["tune"] = {"k": 1}
    _v("tune added", m, "tune")
    m = copy.deepcopy(base)
    m["outputs"] = {"sar": {"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]}}
    _v("outputs.sar added", m, "outputs.sar")
    m = copy.deepcopy(base)
    m["outputs"] = {"b1_plus": {"points": [[0.0, 0.0, 0.0]]}}
    _v("outputs.b1_plus added", m, "outputs.b1_plus")

    fields, expected = [], []
    for label, m, want in variants:
        with pytest.raises(StudyConfigError) as ei:
            _study_from_mapping(m, f"<{label}>")
        _say(capsys, f"  (X) {label}: field {ei.value.field!r} (expected {want!r}); {ei.value}")
        assert type(ei.value) is StudyConfigError
        fields.append(ei.value.field)
        expected.append(want)
    assert fields == expected, fields
    assert calls["n"] == 0


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    # Build cache copied from tests/package/test_wf16_study_tune_capacitors.py:99-109
    # (item 463's fixture, itself item 373's): the first call builds, later
    # equal-geometry calls on the same communicator return that build.
    real_build = FourPortBirdcageGeometry.build
    cache = {"key": None, "built": None, "hits": 0, "real": 0}

    def _build(self, comm_):
        if cache["built"] is not None and cache["key"] == self and comm_ is comm:
            cache["hits"] += 1
            return cache["built"]
        cache["real"] += 1
        out = real_build(self, comm_)
        cache["key"], cache["built"] = self, out
        return out

    out = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        d_t = LOGS / "f_small_hole_pec_tune_s13"
        _fresh_dir(d_t, comm)
        t0 = time.perf_counter()
        run_t = run_study_file(TUNE_FILE, output_dir=d_t, comm=comm)
        out["t_run_t"] = time.perf_counter() - t0
        res_t = run_t.study.result
        t2 = _tune_summary(run_t.study, {"k": 2}, comm)
        assert run_t.study.result is res_t  # no re-solve
        d_c = LOGS / "f_small_hole_pec_capacitors"
        _fresh_dir(d_c, comm)
        t0 = time.perf_counter()
        run_c = run_study_file(CAP_FILE, output_dir=d_c, comm=comm)
        out["t_run_c"] = time.perf_counter() - t0
        comm.Barrier()
        with open(d_c / "summary.json") as fh:
            disk = json.load(fh)
        py = _py_caps()
        t0 = time.perf_counter()
        py_res = py.solve(comm)
        out["t_py"] = time.perf_counter() - t0
    _p0(f"\n[WF-1 s13] build cache: {cache['real']} real build, {cache['hits']} hits (3 Studies)")
    _p0(
        f"[WF-1 s13] tune file: cells {run_t.result.cells}, mesh {run_t.result.mesh_time:.1f} s, "
        f"sweep {run_t.result.sweep_time:.1f} s, run_study_file {out['t_run_t']:.1f} s"
    )
    _p0(
        f"[WF-1 s13] capacitor file: cells {run_c.result.cells}, mesh {run_c.result.mesh_time:.1f} s, "
        f"sweep {run_c.result.sweep_time:.1f} s, run_study_file {out['t_run_c']:.1f} s"
    )
    _p0(f"[WF-1 s13] Python Study: sweep {py_res.sweep_time:.1f} s, total {out['t_py']:.1f} s")
    _p0(
        f"[WF-1 s13] tune k=1 {run_t.summary['tune']['c_mode_farad']:.9e} F, "
        f"k=2 {t2['c_mode_farad']:.9e} F"
    )
    assert cache["real"] == 1 and cache["hits"] == 2
    out.update(run_t=run_t, t2=t2, run_c=run_c, disk=disk, py=py, py_res=py_res)
    return out


@complex_only
def test_wf1_s13_r0_fixture_tie(runs):
    """(R0), (a), (a'): the PEC hole fixture and the transcribed ``C_mode``'s custody."""
    run_t, run_c = runs["run_t"], runs["run_c"]
    assert run_t.summary["cells"] == run_c.summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    means = _class_means(np.asarray(run_t.result.s_matrix, dtype=np.complex128))
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]
    rel = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    c1 = run_t.summary["tune"]["c_mode_farad"]
    a = c1 / TH17_STEP3_C_MODE_RECORD_F - 1.0
    a_prime = ITEM429_C_MODE_F / c1 - 1.0
    _p0(
        f"\n[WF-1 s13 R0] cells {run_t.summary['cells']} / {run_c.summary['cells']}; class means "
        "vs 'pec' record: "
        + ", ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e}; prior 5.052e-11 / 3.878e-11 / 5.875e-11)"
    )
    _p0(
        f"[WF-1 s13 a] C_mode {c1:.9e} F vs TH-17 record {a:+.6e} (band {REDUCTION_BAND:.0e}; "
        "item 442 -2.976552e-4)"
    )
    _p0(
        f"[WF-1 s13 a'] ITEM429_C_MODE_F {ITEM429_C_MODE_F.hex()} vs in-run {float(c1).hex()}: "
        f"{a_prime:+.3e} (band {REDUCTION_BAND:.0e}; item 463 {ITEM463_A_PRIME:+.3e}; "
        f"vs EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.0e}: "
        f"{'inside' if abs(a_prime) <= EXACT_IDENTITY_RTOL else 'outside'}, printed)"
    )
    assert max(rel.values()) <= LEG_D0_REPRODUCTION_BAND
    assert abs(a) <= REDUCTION_BAND
    assert abs(a_prime) <= REDUCTION_BAND


@complex_only
def test_wf1_s13_summary_equals_same_study(runs):
    """(I=): the summary is the same Study's result, bitwise; (c0) ``kept_ports``."""
    run_c, disk = runs["run_c"], runs["disk"]
    s = run_c.summary
    res = run_c.study.result
    exc = run_c.study.excitation
    assert res is run_c.result
    assert len(s["s_matrix"]) == 1 and len(s["s_matrix"][0]) == 1
    for i, row in enumerate(s["s_matrix"]):
        for j, (re, im) in enumerate(row):
            assert (float(re).hex(), float(im).hex()) == _hex(res.s_matrix[i, j])
    for key in ("s_matrix", "kept_ports", "capacitors", "cells"):
        assert json.dumps(disk[key]) == json.dumps(s[key]), key
    for i, row in enumerate(disk["s_matrix"]):
        for j, (re, im) in enumerate(row):
            assert (float(re).hex(), float(im).hex()) == _hex(
                s["s_matrix"][i][j][0] + 1j * s["s_matrix"][i][j][1]
            )
    zc = exc.impedance_ohm(FREQUENCY_64_HZ)
    caps = s["capacitors"]
    _p0(
        f"\n[WF-1 s13 I=] S11 {complex(res.s_matrix[0, 0]):.9e}; kept {s['kept_ports']}; "
        f"capacitors {caps}; keys {sorted(s)}"
    )
    assert s["kept_ports"] == ["P1"]
    assert caps["sheets"] == ["P2", "P3", "P4"]
    assert float(caps["c_f"]).hex() == float(exc.c_f).hex()
    assert float(caps["esr_ohm"]).hex() == float(exc.esr_ohm).hex()
    zs = complex(caps["impedance_ohm"][0], caps["impedance_ohm"][1])
    assert _hex(zs) == _hex(zc) == _hex(res.capacitor_impedance_ohm)
    assert "touchstone" not in s and "z_matrix" not in s


@complex_only
def test_wf1_s13_yaml_equals_python_study(runs):
    """(I-Py): the YAML Study = the Python ``CapacitorSheets`` Study on one mesh."""
    sy = complex(runs["run_c"].result.s_matrix[0, 0])
    sp = complex(runs["py_res"].s_matrix[0, 0])
    assert _spec_diff(runs["run_c"].study, runs["py"]) == []
    r = abs(sy - sp) / abs(sp)
    _p0(
        f"\n[WF-1 s13 I-Py] YAML {sy:.12e} vs Python {sp:.12e}: rel {r:.3e} (band "
        f"{EXACT_IDENTITY_RTOL:.0e}; predicted 0); hex equal {_hex(sy) == _hex(sp)} (printed)"
    )
    assert r <= EXACT_IDENTITY_RTOL


def _route_a(s0, c_f):
    zc = CapacitorSheets(c_f, 0.0, on=SHEETS).impedance_ohm(FREQUENCY_64_HZ)
    # 0-based indices: P2, P3, P4 are rows / columns 1, 2, 3
    return zc, reduce_terminated_ports(s0, float(REFERENCE_IMPEDANCE_OHM), {1: zc, 2: zc, 3: zc})


@complex_only
def test_wf1_s13_route_a_identity(runs):
    """(I-A): the file's S11 equals the Route-A reduction of the tune file's 4x4."""
    run_t, run_c = runs["run_t"], runs["run_c"]
    s0 = np.asarray(run_t.result.s_matrix, dtype=np.complex128)
    assert s0.shape == (4, 4)
    zc = run_c.study.excitation.impedance_ohm(FREQUENCY_64_HZ)
    pa = reduce_terminated_ports(s0, float(REFERENCE_IMPEDANCE_OHM), {1: zc, 2: zc, 3: zc})
    sy = complex(run_c.result.s_matrix[0, 0])
    r = abs(sy - pa[0, 0]) / abs(pa[0, 0])
    _p0(
        f"\n[WF-1 s13 I-A] S11 {sy:.9e} vs reduction {complex(pa[0, 0]):.9e}: rel {r:.6e} "
        f"({r / REDUCTION_BAND:.4f} x REDUCTION_BAND {REDUCTION_BAND:.0e}); item 463 "
        f"{ITEM463_IA:.6e}, S11 {ITEM463_S11:.9e} (printed)"
    )
    assert r <= REDUCTION_BAND


@complex_only
def test_wf1_s13_k2_control(runs):
    """(c1): the k = 2 capacitance on (I-A)'s comparator misses the band."""
    run_t, run_c, t2 = runs["run_t"], runs["run_c"], runs["t2"]
    assert t2["k"] == 2
    s0 = np.asarray(run_t.result.s_matrix, dtype=np.complex128)
    zc2, pa2 = _route_a(s0, t2["c_mode_farad"])
    zc1 = run_c.study.excitation.impedance_ohm(FREQUENCY_64_HZ)
    sy = complex(run_c.result.s_matrix[0, 0])
    r = abs(sy - pa2[0, 0]) / abs(pa2[0, 0])
    _p0(
        f"\n[WF-1 s13 c1] C_2 {t2['c_mode_farad']:.9e} F, z_cap2 {complex(zc2):.6e} vs z_cap1 "
        f"{complex(zc1):.6e}; rel {r:.4e} ({r / REDUCTION_BAND:.1f} x REDUCTION_BAND); "
        f"item 463 {ITEM463_C2:.4e} (printed)"
    )
    assert zc2 != zc1
    assert r > REDUCTION_BAND
