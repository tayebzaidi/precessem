"""`WF-15` step 12 (§9 item 729, §10 U110): the two ungated rod routes refused by field.

No solve.  A rod ``Study`` with ``CapacitorSheets`` (``TypeError``, YAML field
``implant``) or at ``degree != 1`` (``ValueError``, YAML field ``study.degree``) is
refused before any mesh: (Ref) on every entry point that reaches a solve, (NoB) the
builder patched to raise is never called.  CLI: rc [2, 2], no output directory.

Controls: the same ``Study`` with ``implant_rod=None`` passes the check at degree 2, and
degree 1 with no ``CapacitorSheets`` passes it, each reaching a builder patched to a
sentinel; with the refusal patched out, the rod reaches that sentinel too.

`WF-15` step 13 (§9 item 802) lifted the ``CapacitorSheets`` refusal on ``Study.solve``
on the PEC four-port at degree 1 (and its YAML twin); the pins below moved with it.

Not claimed: that either route is wrong.  They are ungated.  ``estimate`` builds a mesh
and solves nothing, so it stays open, and that is pinned here.
"""

from __future__ import annotations

import pytest
from mpi4py import MPI

import fem_em_solver.cli as cli
import fem_em_solver.study as study_module
from fem_em_solver import CapacitorSheets, FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.study_yaml import StudyConfigError

from tests.mesh.test_birdcage_implant_rod import ROD
from tests.package.test_wf15_study_implant_rod import _load, _mapping, _phantom, _rod_yaml
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ

COMM = MPI.COMM_WORLD
CAPS = CapacitorSheets(c_f=1.0e-11, esr_ohm=0.2, on=("P2",))
CAPS_YAML = {"kind": "capacitor_sheets", "c_f": 1.0e-11, "sheets": ["P2"]}
SWEEP_HZ = (64.0e6, 128.0e6)
STUDY = {"frequency_hz": float(STEP2D_FREQUENCY_HZ)}


class _Reached(Exception):
    """The builder's sentinel: the refusal let the call through."""


def _p0(*args):
    if COMM.rank == 0:
        print(*args, flush=True)


def _study(rod, excitation=None, degree=1):
    geometry = FourPortBirdcageGeometry(implant_rod=ROD if rod else None)
    return Study(geometry, _phantom(), excitation or LumpedPortDrive(),
                 frequency_hz=STEP2D_FREQUENCY_HZ, degree=degree)


@pytest.fixture
def builder(monkeypatch):
    """``build`` patched to the sentinel, counting its calls."""
    calls = []

    def _build(self, comm=MPI.COMM_WORLD):
        calls.append(self)
        raise _Reached

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _build)
    return calls


ENTRIES = {
    "solve": lambda s, tmp: s.solve(COMM),
    "solve palace": lambda s, tmp: s.solve(COMM, backend="palace", case_dir=tmp / "case"),
    "solve_sweep": lambda s, tmp: s.solve_sweep(SWEEP_HZ, COMM),
    "sar": lambda s, tmp: s.sar(centres=[(0.0, 0.0, 0.0)], comm=COMM),
    "b1_plus": lambda s, tmp: s.b1_plus(comm=COMM),
    "b1_plus_at": lambda s, tmp: s.b1_plus_at([(0.0, 0.0, 0.0)], comm=COMM),
}


@pytest.mark.parametrize("name", sorted(ENTRIES))
def test_s12_api_refusals_before_any_mesh(name, builder, tmp_path):
    """(Ref), (NoB): every entry point, both refusals, no build.

    `WF-15` step 13 (§9 item 802) moved the ``CAPS`` pin by the 11-02 weekly's ruling
    ("the refusal lifts on exactly the path a step has gated"): ``Study.solve`` on the PEC
    rod at degree 1 is gated there (``test_wf15_study_rod_caps.py``) and reaches the
    builder here; ``solve palace`` passes the Study's check and is refused by the Palace
    writer's own ``implant_rod`` ValueError, before any mesh; every other entry point keeps
    the TypeError.  ``CAPS`` at degree 2 stays refused on every entry point (the degree
    ValueError on ``solve``, the TypeError elsewhere).
    """
    caps = None
    if name == "solve":
        with pytest.raises(_Reached):
            ENTRIES[name](_study(True, CAPS), tmp_path)
    elif name == "solve palace":
        with pytest.raises(ValueError, match="implant_rod") as caps:
            ENTRIES[name](_study(True, CAPS), tmp_path)
    else:
        with pytest.raises(TypeError, match="CapacitorSheets") as caps:
            ENTRIES[name](_study(True, CAPS), tmp_path)
    # Off ``solve`` the CapacitorSheets TypeError still fires first, as on `main`.
    lifted = name in ("solve", "solve palace")
    with pytest.raises(ValueError if lifted else TypeError,
                       match="degree 2" if lifted else "CapacitorSheets") as caps2:
        ENTRIES[name](_study(True, CAPS, degree=2), tmp_path)
    with pytest.raises(ValueError, match="degree 2") as deg:
        ENTRIES[name](_study(True, degree=2), tmp_path)
    assert caps is None or "implant_rod" in str(caps.value)
    assert "implant_rod" in str(deg.value) and "implant_rod" in str(caps2.value)
    assert len(builder) == (1 if name == "solve" else 0)
    _p0(f"\n[WF-15 s12] (Ref) {name}: {caps.value if caps else 'CAPS reached the builder (s13)'}"
        f" | {caps2.value} | {deg.value}; builds {len(builder)}")


def test_s12_tune_and_the_cache(builder):
    """``tune`` solves through ``solve``; a cached field set does not reopen the route."""
    with pytest.raises(ValueError, match="degree 2"):
        _study(True, degree=2).tune(1, comm=COMM)
    cached = _study(True, degree=2)
    cached.fields_result = object()
    with pytest.raises(ValueError, match="degree 2"):
        cached._solve_port_fields(COMM)
    assert builder == []


CONTROLS = {
    "rod-free, degree 2": (lambda: _study(False, degree=2), ("solve", "solve_sweep", "sar")),
    "rod, degree 1": (lambda: _study(True), ("solve", "solve_sweep", "sar")),
    "rod-free, CapacitorSheets": (lambda: _study(False, CAPS), ("solve",)),
}


@pytest.mark.parametrize("name", sorted(CONTROLS))
def test_s12_controls_pass_the_check(name, builder, tmp_path):
    make, entries = CONTROLS[name]
    for entry in entries:
        with pytest.raises(_Reached):
            ENTRIES[entry](make(), tmp_path)
    assert len(builder) == len(entries)
    _p0(f"\n[WF-15 s12] control {name}: reached the builder on {entries}")


def test_s12_the_refusal_is_what_stops_it(builder, monkeypatch, tmp_path):
    """With the check patched out, both refused studies reach the builder on ``solve``."""
    monkeypatch.setattr(study_module, "_refuse_rod_route", lambda study, method: None)
    for study in (_study(True, degree=2), _study(True, CAPS, degree=2)):
        with pytest.raises(_Reached):
            ENTRIES["solve"](study, tmp_path)
    assert len(builder) == 2


def test_s12_estimate_stays_open(builder):
    """``estimate`` solves nothing: a price, not a route (disclosed in the outcome)."""
    with pytest.raises(_Reached):
        _study(True, degree=2).estimate(COMM)
    assert len(builder) == 1


def _doc(study, excitation=None, **extra):
    doc = _mapping(study={**STUDY, **study}, implant=_rod_yaml(), **extra)
    if excitation:
        doc["excitation"] = excitation
    return doc


YAML_CASES = {
    # `WF-15` step 13 (§9 item 802): ``caps`` at degree 1 on the PEC four-port is admitted
    # (test_wf15_study_rod_caps.py); with a tune it meets the YAML's tune refusal, at
    # degree 2 the degree refusal, and on the copper hole route the ``implant`` refusal.
    "caps, tune": (_doc({"degree": 1}, CAPS_YAML, tune={"k": 1}), "tune"),
    "caps, degree 2": (_doc({"degree": 2}, CAPS_YAML), "study.degree"),
    "caps, hole conductor": (_doc({"degree": 1}, CAPS_YAML, conductor={"kind": "leontovich"},
                                  geometry={"kind": "four_port_birdcage", "preset": "f_small",
                                            "as_hole": True}), "implant"),
    "degree 2": (_doc({"degree": 2}), "study.degree"),
    "degree omitted": (_doc({}), "study.degree"),
    "degree 2, sweep": (_doc({"degree": 2, "sweep_hz": list(SWEEP_HZ)}), "study.degree"),
    "degree 2, tune": (_doc({"degree": 2}, tune={"k": 1}), "study.degree"),
}


@pytest.mark.parametrize("name", sorted(YAML_CASES))
def test_s12_yaml_refusals_by_field(name, builder, tmp_path):
    doc, field = YAML_CASES[name]
    with pytest.raises(StudyConfigError) as err:
        _load(tmp_path, "refused", doc)
    assert err.value.field == field, (name, err.value.field)
    assert builder == []
    _p0(f"\n[WF-15 s12] YAML {name}: {err.value.field}: {err.value.reason}")


def test_s12_yaml_controls_parse(tmp_path):
    assert _load(tmp_path, "rod", _doc({"degree": 1})).degree == 1
    free = _doc({"degree": 2})
    del free["implant"]
    assert _load(tmp_path, "free", free).geometry.implant_rod is None
    caps = _doc({"degree": 1}, CAPS_YAML)
    del caps["implant"]
    assert isinstance(_load(tmp_path, "caps", caps).excitation, CapacitorSheets)


def test_s12_study_default_is_what_the_omitted_refusal_assumes():
    """The YAML reads an omitted degree as not 1 because the default is 2."""
    assert _study(False).degree == 1
    assert Study(FourPortBirdcageGeometry(), _phantom(), LumpedPortDrive(),
                 frequency_hz=STEP2D_FREQUENCY_HZ).degree == 2


@pytest.mark.parametrize("name", ["caps, degree 2", "degree 2"])
def test_s12_cli_rc(name, builder, tmp_path):
    import yaml

    path = tmp_path / f"refused-{COMM.rank}.yaml"
    path.write_text(yaml.safe_dump(YAML_CASES[name][0], sort_keys=False))
    out = tmp_path / "out"
    rc = cli.main(["run", str(path), "--output-dir", str(out)])
    rcs = COMM.allgather(rc)
    made = COMM.allgather(out.exists())
    _p0(f"\n[WF-15 s12] CLI {name}: rc {rcs}, output dir {made}")
    assert rcs == [2] * COMM.size
    assert not any(made)
    assert builder == []
