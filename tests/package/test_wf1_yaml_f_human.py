"""`WF-1` step 7a (§9 item 390a, chain U43 part 1): the F-human schema, no mesh.

Schema v0 gains ``geometry.kind: f_human_birdcage`` (hole route only:
``as_hole`` required and ``true``), ``excitation.kind: two_port_quadrature``
(``TwoPortQuadratureDrive``; it needs a ``conductor:`` section) and the
F-human branch of ``run_study_file``, so a YAML file drives `WF-16` step 9's
tuned copper F-human ``Study`` (brief: docs/planning/chunks/WF-1.md, "§9 item
390a").  Nothing here builds a mesh or solves.  Asserted:

(K)    ``load_study`` on the tracked file equals item 382's Python ``Study``
       field by field by ``float.hex``; its generator kwargs equal item 382's
       (K0) comparand; ``cli.main(["check", ...])`` returns 0 on every rank;
(K-c)  the literal SAR lattice (18 lattice points + the c* mirror image) and
       the two ``b1_plus`` points parse bit for bit;
(X)    three tracked malformed files and the in-process one-key mutations
       raise ``StudyConfigError`` naming the field (run / check return 2 on
       every rank, nothing written);
(Run-0) ``run_study_file``'s F-human branch on a counting fake of the tuned
       solve: called once, summary values equal the fake's exactly, no
       Touchstone, the JSON equals the in-memory summary.
Negative controls (asserted, deterministic): (c1) a dropped ``esr_ohm`` diffs
exactly ``['excitation.esr_ohm']``; (c2) a one-ULP ``c_in_f`` diffs exactly
``['excitation.c_in_f']``.

Every helper and constant is imported; the two record constants below are
transcribed with their log lines (records, never bands).
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FHumanBirdcageGeometry,
    Study,
    StudyConfigError,
    TwoPortQuadratureDrive,
    cli,
    load_study,
)
from fem_em_solver import study as study_mod
from fem_em_solver.study_yaml import _load, _study_from_mapping, run_study_file
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _p0,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf16_study_tuned_f_human import _c_f, _study
from tests.package.test_wf16_study_tuned_f_human_b1_sar import H9P_CSTAR_LABEL
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy
from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs
from tests.validation.test_wf10_f_human_copper_esr_root import ESR_OHM
from tests.validation.test_wf10_f_human_tuned_sar10g import AZ17_RECORD_DEG, _lattice

FH_FILE = WF1 / "f_human_tuned_copper.yaml"
MALFORMED_FH = {
    "malformed_f_human_without_hole.yaml": "geometry.as_hole",
    "malformed_quadrature_on_four_port.yaml": "excitation.kind",
    "malformed_quadrature_without_conductor.yaml": "conductor",
}
#: 382b's z_c, printed ``0.000000e+00`` (``.6e`` prints that only for ±0.0):
#: ``20261011T214947Z_WF-16-s9b-W1.log:6383``.  A record, not a band.
Z_C_RECORD_M = 0.0
#: 382's ``az_mirror 56.250000000000 vs JSON 56.250000000000 (|d| 0.000e+00)``:
#: ``20261011T213654Z_WF-16-s9a-W1.log:6369``.  A record, not a band.
AZ_MIRROR_RECORD_DEG = 56.25


def _yf(x) -> str:
    """A float as YAML 1.1 reads it back bit for bit (needs a dot)."""
    r = repr(float(x))
    if "." not in r and "inf" not in r and "nan" not in r:
        r = "%.17e" % float(x)
    return r


def _cstar() -> tuple:
    rows = _lattice(AZ17_RECORD_DEG, Z_C_RECORD_M)
    hits = [row[4] for row in rows if row[0] == H9P_CSTAR_LABEL]
    assert len(hits) == 1, (H9P_CSTAR_LABEL, [row[0] for row in rows])
    return tuple(float(v) for v in hits[0])


def _expected_centres() -> list:
    """The 19 SAR centres: ``_lattice``'s 18 rows in order, then M(c*)."""
    rows = _lattice(AZ17_RECORD_DEG, Z_C_RECORD_M)
    pts = [tuple(float(v) for v in row[4]) for row in rows]
    m = _mirror_xy(np.array([_cstar()], dtype=float), AZ_MIRROR_RECORD_DEG)[0]
    return pts + [tuple(float(v) for v in m)]


def _expected_points() -> list:
    return [(0.0, 0.0, 0.0), _cstar()]


def _literal_lines() -> str:
    """The YAML text of ``c_in_f``, the 19 centres and the two points."""
    lines = [f"  c_in_f: {_yf(_c_f())}", "    centres:"]
    for c in _expected_centres():
        lines.append(f"      - [{', '.join(_yf(v) for v in c)}]")
    lines.append("    points:")
    for q in _expected_points():
        lines.append(f"      - [{', '.join(_yf(v) for v in q)}]")
    return "\n".join(lines)


def _hex3(t):
    return tuple(float(v).hex() for v in t)


def _raises(mapping, source="<mutation>"):
    with pytest.raises(StudyConfigError) as ei:
        _study_from_mapping(mapping, source)
    return ei.value


def test_wf1_f_human_yaml_parse_and_refusals():
    """(K), (K-c), (X), (c1), (c2) — no mesh."""
    import yaml

    comm = MPI.COMM_WORLD
    # ---- (K) bitwise up to the solve --------------------------------------
    s_yaml, s_py = load_study(FH_FILE), _study()
    diff = _spec_diff(s_yaml, s_py)
    _p0(f"\n[WF-1 s7a] (K) diff {diff}")
    _p0(
        f"  c_in_f {s_yaml.excitation.c_in_f.hex()} vs _c_f() {_c_f().hex()}; "
        f"esr_ohm {s_yaml.excitation.esr_ohm.hex()} vs ESR_OHM {float(ESR_OHM).hex()}"
    )
    _p0(
        f"  excitation {s_yaml.excitation!r}; degree {s_yaml.degree}; conductor {s_yaml.conductor!r}"
    )
    assert type(s_yaml.geometry) is FHumanBirdcageGeometry
    assert type(s_yaml.excitation) is TwoPortQuadratureDrive
    assert isinstance(s_yaml, Study)
    assert diff == [], diff
    assert s_yaml.excitation.c_in_f.hex() == _c_f().hex()
    assert s_yaml.excitation.esr_ohm.hex() == float(ESR_OHM).hex()
    kw_yaml, kw_py = s_yaml.geometry.generator_kwargs(), _f_human_kwargs(as_hole=True)
    _p0(f"  (K0) generator_kwargs equal: {kw_yaml == kw_py}")
    assert kw_yaml == kw_py
    rc = cli.main(["check", str(FH_FILE)])
    _p0(f"  (K) cli check rc (min, max) {_rc_all(rc, comm)}")
    assert _rc_all(rc, comm) == (0, 0)

    # ---- (K-c) the literal lattice -----------------------------------------
    _, kind, outputs, cooling = _load(FH_FILE)
    assert kind == "f_human_birdcage" and cooling is None
    exp_c, got_c = _expected_centres(), outputs["sar"]["centres"]
    exp_p, got_p = _expected_points(), outputs["b1_plus"]["points"]
    for i, e in enumerate(exp_c):
        _p0(f"  (K-c) centre[{i:2d}] expected {[repr(v) for v in e]}")
    n_eq = (
        sum(1 for e, g in zip(exp_c, got_c) if _hex3(e) == _hex3(g))
        if len(exp_c) == len(got_c)
        else -1
    )
    _p0(
        f"  (K-c) centres equal by hex: {n_eq} of {len(exp_c)} (parsed {len(got_c)}); "
        f"mass_g {outputs['sar']['mass_g'].hex()}; points {[_hex3(q) for q in got_p]}"
    )
    assert len(exp_c) == 19 and len(got_c) == 19
    assert n_eq == 19
    assert outputs["sar"]["mass_g"].hex() == (10.0).hex()
    assert len(got_p) == 2 and [_hex3(q) for q in got_p] == [_hex3(q) for q in exp_p]

    # ---- (X) tracked malformed files ----------------------------------------
    out = LOGS / "f_human_malformed"
    for name, expected in MALFORMED_FH.items():
        path = WF1 / name
        with pytest.raises(StudyConfigError) as ei:
            load_study(path)
        exc = ei.value
        _fresh_dir(out, comm)
        run_all = _rc_all(cli.main(["run", str(path), "--output-dir", str(out)]), comm)
        check_all = _rc_all(cli.main(["check", str(path)]), comm)
        exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
        _p0(
            f"  (X) {name}: field {exc.field!r} (expected {expected!r}); run rc {run_all}, "
            f"check rc {check_all}, out exists {exists}; {exc}"
        )
        assert type(exc) is StudyConfigError
        assert exc.field == expected
        assert run_all == (2, 2) and check_all == (2, 2)
        assert not exists

    # ---- (X) in-process one-key mutations -----------------------------------
    base = _mapping(FH_FILE)

    def mut(fn):
        m = copy.deepcopy(base)
        fn(m)
        return m

    def _del(sec, key):
        def f(m):
            del m[sec][key]

        return f

    def _set(sec, key, value):
        def f(m):
            m[sec][key] = value

        return f

    def _set_out(key, value):
        def f(m):
            m["outputs"][key] = value

        return f

    def _set_sar_centres(m):
        m["outputs"]["sar"]["centres"] = {"ring_radius_m": 0.045}

    def _add_thermal(m):
        m["thermal"] = {"k_w_per_m_k": 0.6, "h_w_per_m2_k": 10.0, "t_inf_c": 20.0}

    cases = [
        ("unknown geometry key", _set("geometry", "ring_raduis", 0.15), "geometry.ring_raduis"),
        (
            "emit_leg_ports true",
            _set("geometry", "emit_leg_ports", True),
            "geometry.emit_leg_ports",
        ),
        ("as_hole false", _set("geometry", "as_hole", False), "geometry.as_hole"),
        ("as_hole 1", _set("geometry", "as_hole", 1), "geometry.as_hole"),
        ("lumped_port_drive", _set("excitation", "kind", "lumped_port_drive"), "excitation.kind"),
        ("unknown excitation key", _set("excitation", "c_in_pf", 118.2437), "excitation.c_in_pf"),
        ("c_in_f removed", _del("excitation", "c_in_f"), "excitation.c_in_f"),
        ("c_in_f 0.0", _set("excitation", "c_in_f", 0.0), "excitation.c_in_f"),
        ("c_in_f string", _set("excitation", "c_in_f", "1.18e-10"), "excitation.c_in_f"),
        ("esr_ohm -0.1", _set("excitation", "esr_ohm", -0.1), "excitation.esr_ohm"),
        ("feeds [P17]", _set("excitation", "feeds", ["P17"]), "excitation.feeds"),
        ("feeds [P17, P17]", _set("excitation", "feeds", ["P17", "P17"]), "excitation.feeds"),
        ("feeds [P17, P1]", _set("excitation", "feeds", ["P17", "P1"]), "excitation.feeds[1]"),
        ("z0_ohm 0.0", _set("excitation", "z0_ohm", 0.0), "excitation.z0_ohm"),
        ("thermal section", _add_thermal, "thermal"),
        ("ring_radius_m centres", _set_sar_centres, "outputs.sar.centres"),
        (
            "probes",
            _set_out("probes", [{"centre": [0.0, 0.0, 0.0], "radius": 0.01}]),
            "outputs.probes",
        ),
    ]
    for label, fn, expected in cases:
        exc = _raises(mut(fn), f"<{label}>")
        _p0(f"  (X) {label}: field {exc.field!r} (expected {expected!r}); {exc.reason}")
        assert type(exc) is StudyConfigError, label
        assert exc.field == expected, (label, exc.field)
    _p0(f"  (X) {len(cases)} of {len(cases)} mutations named their field")

    # ---- (c1) dropped key / (c2) one ULP ------------------------------------
    m1 = mut(_del("excitation", "esr_ohm"))
    s1 = _study_from_mapping(m1)
    d1 = _spec_diff(s1, s_py)
    default_esr = TwoPortQuadratureDrive(_c_f()).esr_ohm
    _p0(
        f"  (c1) dropped esr_ohm diff {d1}; esr {s1.excitation.esr_ohm.hex()} vs default "
        f"{float(default_esr).hex()}"
    )
    assert d1 == ["excitation.esr_ohm"], d1
    assert s1.excitation.esr_ohm.hex() == float(default_esr).hex()
    ulp = yaml.safe_load(repr(float(np.nextafter(_c_f(), 1.0))))
    assert type(ulp) is float and ulp != _c_f()
    s2 = _study_from_mapping(mut(_set("excitation", "c_in_f", ulp)))
    d2 = _spec_diff(s2, s_py)
    _p0(f"  (c2) one-ULP c_in_f {ulp!r} diff {d2}")
    assert d2 == ["excitation.c_in_f"], d2


def test_wf1_f_human_yaml_run_branch_no_mesh():
    """(Run-0): ``run_study_file``'s F-human branch on a fake tuned result."""
    import yaml

    comm = MPI.COMM_WORLD
    s_node = np.array([[0.125 + 0.25j, -0.375 + 0.5j], [-0.375 + 0.5j, 0.625 - 0.75j]])
    y_in = {
        "ccw": np.array([0.01 + 0.02j, 0.03 - 0.04j]),
        "cw": np.array([-0.05 + 0.06j, 0.07 + 0.08j]),
    }
    gap, feed, vs = complex(0.2, -21.03), complex(7.6, -17.7), complex(1.02, 2.38)
    calls = []

    def fake(study, comm_=MPI.COMM_WORLD):
        calls.append(1)
        study.result = SimpleNamespace(
            cells=369387,
            s_node=s_node,
            y_in=y_in,
            feeds=("P17", "P21"),
            az_feeds=(11.25, 101.25),
            az_mirror=56.25,
            gap_impedance_ohm=gap,
            feed_impedance_ohm=feed,
            source_voltage_v=vs,
        )
        return study.result

    def boom(*args, **kwargs):
        raise AssertionError("no mesh / field solve in item 390a")

    run_dir = LOGS / "f_human_run0"
    _fresh_dir(run_dir, comm)
    yaml_path = run_dir / "f_human_run0.yaml"
    if comm.rank == 0:
        m = _mapping(FH_FILE)
        del m["outputs"]
        run_dir.mkdir(parents=True, exist_ok=True)
        yaml_path.write_text(yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    out = run_dir / "out"
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(study_mod, "_study_solve_tuned_quadrature", fake)
        mp.setattr(study_mod.Study, "_solve_port_fields", boom)
        mp.setattr(study_mod.FHumanBirdcageGeometry, "build", boom)
        run = run_study_file(yaml_path, output_dir=out, comm=comm)
    summary = run.summary
    n_calls = comm.allreduce(len(calls), op=MPI.MAX), comm.allreduce(len(calls), op=MPI.MIN)
    _p0(f"\n[WF-1 s7a] (Run-0) fake calls (max, min) over ranks {n_calls}; keys {sorted(summary)}")
    _p0(f"  summary excitation {summary['excitation']}")
    assert n_calls == (1, 1)

    exp_s = [[[0.125, 0.25], [-0.375, 0.5]], [[-0.375, 0.5], [0.625, -0.75]]]
    exp_y = {"ccw": [[0.01, 0.02], [0.03, -0.04]], "cw": [[-0.05, 0.06], [0.07, 0.08]]}
    exp_exc = {
        "kind": "two_port_quadrature",
        "c_in_f": _c_f(),
        "esr_ohm": float(ESR_OHM),
        "feeds": ["P17", "P21"],
        "z0_ohm": 50.0,
        "gap_impedance_ohm": [0.2, -21.03],
        "feed_impedance_ohm": [7.6, -17.7],
        "source_voltage_v": [1.02, 2.38],
    }
    assert summary["family"] == "f_human_birdcage"
    assert summary["cells"] == 369387
    assert summary["s_node"] == exp_s
    assert summary["y_in"] == exp_y
    assert summary["az_feeds"] == [11.25, 101.25]
    assert summary["az_mirror"] == 56.25
    assert summary["excitation"] == exp_exc
    assert summary["excitation"]["c_in_f"].hex() == _c_f().hex()
    assert summary["degree"] == 1 and summary["ranks"] == comm.size
    assert summary["conductor"] == {"kind": "leontovich", "sigma_s_per_m": 58000000.0}
    assert "touchstone" not in summary and "z_matrix" not in summary
    assert "touchstone" not in run.paths
    assert "sar" not in summary and "b1_plus" not in summary
    same = True
    if comm.rank == 0:
        with open(out / "summary.json") as fh:
            same = json.load(fh) == summary
    same = comm.allreduce(bool(same), op=MPI.LAND)
    _p0(f"  (Run-0) written JSON == in-memory summary: {same}")
    assert same
    assert Path(run.paths["summary"]).name == "summary.json"
