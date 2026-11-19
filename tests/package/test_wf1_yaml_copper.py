"""`WF-1` step 4 (§9 item 343, chain U32): ``conductor:`` + ``geometry.as_hole``.

Schema v0 gains an optional top-level ``conductor:`` section (``kind:
leontovich``, optional ``sigma_s_per_m``) and the one settable ``bool``
geometry field ``as_hole``, so ``run_study_file`` drives `WF-16` step 7's copper
port Study (§9 item 324) at degree 1, 64 MHz, ``-n 2`` (brief:
docs/planning/chunks/WF-1.md, "§9 item 343").  Asserted:

(K)  ``load_study("f_small_copper.yaml")`` equals the Python-built copper Study
     field by field by ``float.hex`` (``_spec_diff`` extended to ``conductor``);
     the PEC control's diff is exactly ``['conductor']``; an omitted
     ``sigma_s_per_m`` parses to ``LeontovichWall()``'s default;
(R)  the YAML Cu class means reproduce ``ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]
     ["cu"]`` at ``LEG_D0_REPRODUCTION_BAND``; cells equal
     ``ANS6_EXACT_STRIP_HOLE_CELL_RECORD``; the JSON ``s_matrix`` equals the
     result's exactly;
(b)  the copper power identity on the YAML-solved result: ``p_surf > 0`` and
     residual <= ``DISCRETE_IDENTITY_RTOL``;
(g)  `TH-14`'s gates on the YAML Cu and PEC-control sweeps;
(I)  an independently solved Python Cu Study within ``LEG_D0_REPRODUCTION_BAND``
     of the YAML run (rule (j): bitwise stops at the solve);
(X)  three malformed files raise ``StudyConfigError`` naming ``conductor`` /
     ``conductor`` / ``conductor.sigma_s_per_m``; ``cli.main`` run / check
     return 2 on every rank, nothing written.
Negative control (asserted): (c1) the PEC-control file reproduces the record's
``"pec"`` column and misses the Cu record by > ``NEGATIVE_CONTROL_FLOOR``.

Every band, record and helper is imported, none restated.
"""

from __future__ import annotations

import copy
import json
import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import (
    LeontovichWall,
    LumpedPortDrive,
    Study,
    StudyConfigError,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.study_yaml import _study_from_mapping
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf16_study_copper import (
    NEGATIVE_CONTROL_FLOOR,
    _geometry,
    _identity,
    _rel_means,
)
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KEY = "64"
COPPER_FILE = WF1 / "f_small_copper.yaml"
PEC_CONTROL_FILE = WF1 / "f_small_hole_pec_control.yaml"
MALFORMED_CONDUCTOR = {
    "malformed_conductor_on_sphere.yaml": "conductor",
    "malformed_conductor_without_hole.yaml": "conductor",
    "malformed_conductor_sigma_zero.yaml": "conductor.sigma_s_per_m",
}


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _py_copper(conductor="cu"):
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    cond = LeontovichWall(sigma_s_per_m=study_mod.COPPER_SIGMA) if conductor == "cu" else None
    return Study(
        _geometry(),
        phantom,
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        conductor=cond,
    )


def _fmt(means):
    return ", ".join(f"{c} {means[c].real:+.9e}{means[c].imag:+.9e}j" for c in CLASSES)


def test_wf1_copper_yaml_parse_and_refusals():
    """(K) and (X): no mesh, no solve."""
    comm = MPI.COMM_WORLD
    s_yaml, s_py = load_study(COPPER_FILE), _py_copper()
    d_cu = _spec_diff(s_yaml, s_py)
    d_ctl = _spec_diff(load_study(PEC_CONTROL_FILE), s_py)
    g = s_yaml.geometry
    _p0(
        f"\n[WF-1 s4] (K) copper diff {d_cu}; PEC-control diff {d_ctl}; conductor "
        f"{s_yaml.conductor!r} (sigma hex {s_yaml.conductor.sigma_s_per_m.hex()} vs "
        f"COPPER_SIGMA {study_mod.COPPER_SIGMA.hex()}); as_hole {g.as_hole!r}; "
        f"strip_fraction {g.strip_fraction!r} (GATED_WIDTH_FRACTION "
        f"{study_mod.GATED_WIDTH_FRACTION!r})"
    )
    assert d_cu == [], d_cu
    assert d_ctl == ["conductor"], d_ctl
    assert g.as_hole is True
    assert g.strip_fraction == study_mod.GATED_WIDTH_FRACTION

    # Omitted sigma_s_per_m: LeontovichWall's own default, never a loader one.
    m = copy.deepcopy(_mapping(COPPER_FILE))
    del m["conductor"]["sigma_s_per_m"]
    c_def = _study_from_mapping(m, "<omitted>").conductor
    _p0(f"  (K) omitted sigma_s_per_m -> {c_def!r} (default {LeontovichWall()!r})")
    assert type(c_def) is LeontovichWall
    assert c_def.sigma_s_per_m.hex() == LeontovichWall().sigma_s_per_m.hex()

    # as_hole must be a YAML bool: an int is refused, not coerced.
    m_int = copy.deepcopy(_mapping(COPPER_FILE))
    m_int["geometry"]["as_hole"] = 1
    with pytest.raises(StudyConfigError) as ei:
        _study_from_mapping(m_int, "<as_hole int>")
    _p0(f"  (K) as_hole: 1 -> field {ei.value.field!r}; {ei.value}")
    assert ei.value.field == "geometry.as_hole"

    out = LOGS / "malformed_conductor"
    fields = []
    for name, expected in MALFORMED_CONDUCTOR.items():
        path = WF1 / name
        with pytest.raises(StudyConfigError) as ei:
            load_study(path)
        exc = ei.value
        _p0(f"  (X) {name}: field {exc.field!r} (expected {expected!r}); {exc}")
        assert type(exc) is StudyConfigError
        assert exc.field == expected
        assert exc.source.endswith(name)
        fields.append(exc.field)
        _fresh_dir(out, comm)
        rc_run = cli.main(["run", str(path), "--output-dir", str(out)])
        rc_check = cli.main(["check", str(path)])
        run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
        exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
        _p0(f"  (X) {name}: run rc (min, max) {run_all}, check rc {check_all}, out exists {exists}")
        assert run_all == (2, 2)
        assert check_all == (2, 2)
        assert not exists
    assert fields == list(MALFORMED_CONDUCTOR.values()), fields


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"(R) is asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    assert th14.FREQUENCIES_HZ[KEY] == FREQUENCY_64_HZ
    out = {}
    for col, path in (("cu", COPPER_FILE), ("pec", PEC_CONTROL_FILE)):
        d = LOGS / path.stem
        _fresh_dir(d, comm)
        t0 = time.perf_counter()
        run = run_study_file(path, output_dir=d, comm=comm)
        t = time.perf_counter() - t0
        with open(d / "summary.json") as fh:
            summary = json.load(fh)
        means = _class_means(run.result.s_matrix)
        out[col] = {
            "run": run,
            "summary": summary,
            "means": means,
            "t": t,
            "rec": th14._sweep_record(run.result.sweep),
        }
        _p0(
            f"\n[WF-1 s4] YAML {col} ({path.name}): cells {run.result.cells}, run_study_file "
            f"{t:.1f} s, summary conductor {summary.get('conductor', '<absent>')}"
        )
        _p0(f"    {col} class means (10 digits): {_fmt(means)}")
    t0 = time.perf_counter()
    py = _py_copper()
    py_res = py.solve(comm)
    out["py"] = {"means": _class_means(py_res.s_matrix), "t": time.perf_counter() - t0}
    _p0(f"[WF-1 s4] Python Cu Study {out['py']['t']:.1f} s; class means {_fmt(out['py']['means'])}")
    t0 = time.perf_counter()
    out["power"] = _identity(out["cu"]["run"].result, out["cu"]["run"].study.omega, comm)
    _p0(f"[WF-1 s4] identity field solve {time.perf_counter() - t0:.1f} s; power {out['power']}")
    for col in ("cu", "pec"):
        rel = _rel_means(out[col]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY][col])
        _p0(
            f"[WF-1 s4] (R) YAML {col} vs record {col!r}, relative: "
            + "  ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
            + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
        )
    rel_c1 = _rel_means(out["pec"]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    _p0(
        "[WF-1 s4] (c1) PEC control vs Cu record, relative (PREDICTED 6.420e-4 / 2.811e-4 / "
        f"4.861e-4; floor {NEGATIVE_CONTROL_FLOOR:.0e}): "
        + "  ".join(f"{c} {rel_c1[c]:.3e}" for c in CLASSES)
    )
    rel_i = _rel_means(out["py"]["means"], out["cu"]["means"])
    _p0(
        "[WF-1 s4] (I) Python Cu vs YAML Cu, relative: "
        + "  ".join(f"{c} {rel_i[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    for col in ("cu", "pec"):
        sw = out[col]["rec"]
        _p0(
            f"[WF-1 s4] (g) {col}: reciprocity {sw['reciprocity']:.3e}, max sigma "
            f"{float(np.max(sw['sigma'])):.12f}, spreads {sw['spreads']}"
        )
    return out


@complex_only
def test_wf1_copper_yaml_reproduces_record(runs):
    """(R): class means, cell record, JSON summary."""
    cu = runs["cu"]
    run, summary = cu["run"], cu["summary"]
    rel = _rel_means(cu["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])
    assert run.result.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    s = run.result.s_matrix
    want = [[[float(z.real), float(z.imag)] for z in row] for row in s]
    assert summary["s_matrix"] == want
    assert run.summary["s_matrix"] == want
    assert summary["conductor"] == {"kind": "leontovich", "sigma_s_per_m": study_mod.COPPER_SIGMA}
    assert summary["family"] == "four_port_birdcage" and summary["degree"] == 1
    assert "conductor" not in runs["pec"]["summary"]


@complex_only
def test_wf1_copper_yaml_power_identity(runs):
    """(b) on the YAML-solved result."""
    p = runs["power"]
    assert p["p_surf"] > 0.0, p
    assert p["residual"] <= DISCRETE_IDENTITY_RTOL, p


@complex_only
@pytest.mark.parametrize("col", ["cu", "pec"])
def test_wf1_copper_yaml_th14_gates(runs, col):
    """(g)"""
    sw = runs[col]["rec"]
    assert sw["reciprocity"] <= RECIPROCITY_BAND, (col, sw["reciprocity"])
    assert float(np.max(sw["sigma"])) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (col, sw["sigma"])
    for name, value in sw["spreads"].items():
        assert value <= ADJACENT_SPREAD_BAND, (col, name, value)


@complex_only
def test_wf1_copper_yaml_vs_python_study(runs):
    """(I) at the imported band (rule (j))."""
    rel = _rel_means(runs["py"]["means"], runs["cu"]["means"])
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])


@complex_only
def test_wf1_copper_yaml_negative_control(runs):
    """(c1): conductor dropped reproduces "pec" and misses "cu"."""
    rel_pec = _rel_means(runs["pec"]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["pec"])
    rel_cu = _rel_means(runs["pec"]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    for c in CLASSES:
        assert rel_pec[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel_pec[c])
        assert rel_cu[c] > NEGATIVE_CONTROL_FLOOR, (c, rel_cu[c])
