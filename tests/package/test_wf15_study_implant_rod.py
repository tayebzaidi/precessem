"""`WF-15` step 10 (§9 item 709, §10 U97′): the rod in the coil-driven F-small phantom
through the ``Study`` — ``FourPortBirdcageGeometry(implant_rod=ROD)``, degree 1, 64 MHz,
the four-port ``LumpedPortDrive``, on item 235's exact-rectangle rod record mesh.

Knob ``FEM_EM_WF15_S10`` = ``1``; unset, the coding tests run (no mesh) and the solve
node skips.  Every band is imported; the one constant set here is the cell record.

Asserted: the ``Study``'s cells, sheets and S equal the direct
``_build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`` + ``build_four_port_sweep``
at ``IDENTITY_BAND`` in one window; (i) ``RECIPROCITY_BAND``; (ii) σ_max ≤ 1 +
``PASSIVITY_SIGMA_TOLERANCE`` and column power ≤ 1; (M) σ_y pairs ≤
``ADJACENT_SPREAD_BAND``; (E) rod wall ≤ ``PEC_WALL_MAX_E``.  Controls: the off ``Study``
reproduces ``ANS4_EXACT_STRIP_CELL_RECORD``; tag-402 wall count > 0 with the rod, 0 off.
Printed, never asserted: (PP), never gated on this mesh; the tip 1 g Q per single-port
drive (``Study.sar`` has no ccw, and no record exists on this mesh).

Not claimed: a two-port drive on F-small, the 1 g rung, a tip ratio, (PP), a published
comparand, a kept rod, a clinical rod.
"""

from __future__ import annotations

import dataclasses
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_module
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.study_yaml import StudyConfigError, load_study

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD, _build_rod
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    STEP2D_FREQUENCY_HZ,
    _exact_shares,
)
from tests.validation.test_imp3_palace_twin import _sweep
from tests.validation.test_imp3_rod_birdcage import (
    PEC_WALL_MAX_E,
    PORT_IDS,
    TIPS,
    _matrix_quantities,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_imp3_rod_heating_phantom import _parent_wall_count
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KNOB = "FEM_EM_WF15_S10"
#: Item 235's exact-rectangle rod record mesh, the 0.10 image at `-n 2`: "arm rod:
#: 168507 cells" (`20261005T035322Z_IMP-3-item235-W6-ours-p1.log.gz:3731`).  Not
#: `ROD_CELL_RECORD_W1` (ragged) or `KEPT_ROD_CELL_RECORD` (kept).
ROD_EXACT_CELL_RECORD = 168_507
RECORD_WIDTH = 2

needs_knob = pytest.mark.skipif(os.environ.get(KNOB) != "1", reason=f"{KNOB} unset")
coding = pytest.mark.skipif(os.environ.get(KNOB) == "1", reason=f"{KNOB} set: coding tests skip")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, STEP2D_FREQUENCY_HZ)


def _study(geometry):
    return Study(geometry, _phantom(), LumpedPortDrive(), frequency_hz=STEP2D_FREQUENCY_HZ, degree=1)


class _Stop(Exception):
    pass


def _kwargs_of(monkeypatch, geometry):
    seen = {}

    def _record(**kw):
        seen.update(kw)
        raise _Stop

    monkeypatch.setattr(study_module.MeshGenerator, "birdcage_port_domain", staticmethod(_record))
    with pytest.raises(_Stop):
        geometry.build(MPI.COMM_WORLD)
    return seen


def _same_kwargs(a, b):
    assert set(a) == set(b), set(a) ^ set(b)
    for key in a:
        if isinstance(a[key], np.ndarray):
            assert np.array_equal(a[key], b[key]), key
        else:
            assert a[key] == b[key], key


@coding
def test_s10_field_equality_and_hash():
    """The default equals ``implant_rod=None``; a dict field stays out of the hash."""
    default = FourPortBirdcageGeometry()
    assert default == FourPortBirdcageGeometry(implant_rod=None)
    assert default == FourPortBirdcageGeometry.f_small()
    assert hash(default) == hash(FourPortBirdcageGeometry(implant_rod=None))
    rod = FourPortBirdcageGeometry(implant_rod=ROD)
    assert rod != default
    assert rod == FourPortBirdcageGeometry(implant_rod=dict(ROD))
    assert isinstance(hash(rod), int)
    assert dataclasses.replace(rod, implant_rod=None) == default


@coding
def test_s10_generator_kwargs_identity(monkeypatch):
    """`None` forwards nothing: the generator call is unchanged, with no mesh; the rod
    adds exactly one key; κ is reached by neither."""
    base = _kwargs_of(monkeypatch, FourPortBirdcageGeometry())
    none = _kwargs_of(monkeypatch, FourPortBirdcageGeometry(implant_rod=None))
    rod = _kwargs_of(monkeypatch, FourPortBirdcageGeometry(implant_rod=ROD))
    _same_kwargs(base, none)
    assert "implant_rod" not in base
    assert set(rod) - set(base) == {"implant_rod"} and rod["implant_rod"] == ROD
    assert rod["implant_rod"] is not ROD
    _same_kwargs({k: v for k, v in rod.items() if k != "implant_rod"}, base)
    assert not any("kappa" in k for k in rod), sorted(rod)
    _p0(f"\n[WF-15 s10] kwargs: {len(base)} keys unchanged by implant_rod=None; the rod adds one")


def _mapping(**extra):
    out = {
        "femem": 1,
        "study": {"frequency_hz": float(STEP2D_FREQUENCY_HZ), "degree": 1},
        "geometry": {"kind": "four_port_birdcage", "preset": "f_small"},
        "material": {"sigma": float(SALINE_SIGMA), "epsilon_r": float(SALINE_EPSILON_R)},
        "excitation": {"kind": "lumped_port_drive"},
    }
    out.update(extra)
    return out


def _rod_yaml():
    return {"rod": {"radius_m": ROD["radius"], "length_m": ROD["length"],
                    "centre_m": list(ROD["centre"]), "axis": list(ROD["axis"])}}


def _load(tmp_path, name, mapping):
    import yaml

    comm = MPI.COMM_WORLD
    path = tmp_path / f"{name}-{comm.rank}.yaml"
    path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    return load_study(path)


@coding
def test_s10_yaml_section_and_refusals(tmp_path):
    """`implant:` parses to the Python mapping exactly; it is refused by field on every
    other kind, beside `load_result:`, and on a malformed section."""
    st = _load(tmp_path, "good", _mapping(implant=_rod_yaml()))
    assert st.geometry == FourPortBirdcageGeometry(implant_rod=ROD)
    for key in ("radius", "length"):
        assert float(st.geometry.implant_rod[key]).hex() == float(ROD[key]).hex(), key
    for key in ("centre", "axis"):
        assert st.geometry.implant_rod[key] == tuple(float(v) for v in ROD[key]), key
    assert _load(tmp_path, "none", _mapping()).geometry.implant_rod is None

    sphere = {"kind": "sphere_in_box", "sphere_radius": 0.05, "box_half_width": 0.1,
              "resolution_sphere": 0.02, "resolution_far": 0.04}
    ball = {"kind": "thermal_ball", "radius": 0.05, "resolution": 0.01}
    missing = _rod_yaml()
    del missing["rod"]["axis"]
    bad_triple = _rod_yaml()
    bad_triple["rod"]["centre_m"] = [0.0, 0.0]
    unknown = _rod_yaml()
    unknown["rod"]["keep_cells"] = True
    cases = {
        "sphere": (_mapping(geometry=sphere, implant=_rod_yaml()), "implant"),
        "human": (_mapping(geometry={"kind": "f_human_birdcage", "as_hole": True},
                           implant=_rod_yaml()), "implant"),
        "clinical": (_mapping(geometry={"kind": "f_clinical_birdcage"}, implant=_rod_yaml()), "implant"),
        "ball": ({"femem": 1, "geometry": ball, "excitation": {}, "thermal": {},
                  "implant": _rod_yaml()}, "implant"),
        "loaded": (_mapping(implant=_rod_yaml(), load_result="x.npz", tune={}), "implant"),
        "not_a_mapping": (_mapping(implant="rod"), "implant"),
        "no_rod": (_mapping(implant={}), "implant.rod"),
        "extra": (_mapping(implant={**_rod_yaml(), "spheroid": {}}), "implant.spheroid"),
        "missing": (_mapping(implant=missing), "implant.rod.axis"),
        "triple": (_mapping(implant=bad_triple), "implant.rod.centre_m"),
        "unknown": (_mapping(implant=unknown), "implant.rod.keep_cells"),
        "geometry_key": (_mapping(geometry={"kind": "four_port_birdcage", "implant_rod": {}}),
                         "geometry.implant_rod"),
    }
    for name, (mapping, field) in cases.items():
        with pytest.raises(StudyConfigError) as err:
            _load(tmp_path, name, mapping)
        _p0(f"[WF-15 s10] (Ref) {name}: {err.value.field}: {err.value.reason}")
        assert err.value.field == field, (name, err.value.field)


@coding
def test_s10_palace_refuses_a_rod(tmp_path):
    """The Palace case carries no rod-wall boundary, so it refuses before any mesh."""
    from fem_em_solver.backends.palace import solve_palace

    with pytest.raises(ValueError, match="implant_rod"):
        solve_palace(_study(FourPortBirdcageGeometry(implant_rod=ROD)), MPI.COMM_WORLD,
                     case_dir=tmp_path)


@pytest.fixture(scope="module")
def s10():
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_WIDTH, f"the record holds at -n {RECORD_WIDTH}, got -n {comm.size}"
    assert ROD_EXACT_CELL_RECORD != ANS4_EXACT_STRIP_CELL_RECORD
    captured = {}
    real_build = FourPortBirdcageGeometry.build

    def _capturing_build(self, comm_=MPI.COMM_WORLD):
        out = real_build(self, comm_)
        captured["rod" if self.implant_rod else "off"] = out
        return out

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _capturing_build)
        t0 = time.perf_counter()
        rod_study = _study(FourPortBirdcageGeometry(implant_rod=ROD))
        rod = rod_study._solve_port_fields(comm)
        t_rod = time.perf_counter() - t0
        t0 = time.perf_counter()
        centres = list(TIPS.values())
        rod_sar = rod_study.sar(1.0, centres=centres, comm=comm)
        t_sar = time.perf_counter() - t0

        t0 = time.perf_counter()
        built = _build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
        case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
        direct = _sweep(case, keep=True)
        t_direct = time.perf_counter() - t0
        omega = 2.0 * np.pi * float(STEP2D_FREQUENCY_HZ)
        shares = {pid: _exact_shares(case, {"driven": pid, "omega": omega, "fields": direct.fields[pid]})
                  for pid in PORT_IDS}

        t0 = time.perf_counter()
        off_study = _study(FourPortBirdcageGeometry())
        off = off_study._solve_port_fields(comm)
        off_sar = off_study.sar(1.0, centres=centres, comm=comm)
        t_off = time.perf_counter() - t0
    _p0(f"\n[WF-15 s10] -n {comm.size}: rod Study {t_rod:.1f} s (mesh {rod.mesh_time:.1f}, sweep "
        f"{rod.sweep_time:.1f}), its 1 g reads {t_sar:.1f} s; direct {t_direct:.1f} s; off Study "
        f"{t_off:.1f} s")
    return {"rod": rod, "case": case, "direct": direct, "shares": shares, "off": off,
            "rod_sar": rod_sar, "off_sar": off_sar, "captured": captured, "built": built}


@complex_only
@needs_knob
def test_s10_all_gates(s10):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    comm = MPI.COMM_WORLD
    rod, direct, case = s10["rod"], s10["direct"], s10["case"]
    mq = _matrix_quantities(rod)
    s_st, s_di = np.asarray(rod.s_matrix), np.asarray(direct.s_matrix)
    rel_s = float(np.linalg.norm(s_st - s_di) / np.linalg.norm(s_di))
    gates = {
        "cells_record": rod.cells == ROD_EXACT_CELL_RECORD,
        "cells_direct": rod.cells == case["cells"],
        "S_direct": rel_s <= IDENTITY_BAND,
        "i": mq["recip"] <= RECIPROCITY_BAND,
        "ii_sigma": float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE,
        "ii_colpow": bool(np.all(mq["colpow"] <= 1.0)),
        "M": max(mq["mirror"]) <= ADJACENT_SPREAD_BAND,
    }
    _p0(f"[WF-15 s10] cells Study {rod.cells}, direct {case['cells']}, record {ROD_EXACT_CELL_RECORD}; "
        f"||S_study - S_direct||/||S_direct|| {rel_s:.3e} (IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED)")
    _p0(f"[WF-15 s10] (i) {mq['recip']:.3e} (band {RECIPROCITY_BAND:.0e}); (ii) sigma_max "
        f"{float(mq['sigma'][0]):.12f}, colpow {np.array2string(mq['colpow'], precision=12)}; (M) "
        + ", ".join(f"{v:.3e}" for v in mq["mirror"]) + f" (band {ADJACENT_SPREAD_BAND})")
    assert len(rod.sheets) == len(case["sheets"]) == 4
    for a, b in zip(rod.sheets, case["sheets"]):
        gates[f"sheet {a['tag']} facets"] = a["tag"] == b["tag"] and a["facets"] == b["facets"]
        for key in ("h", "w", "area"):
            gates[f"sheet {a['tag']} {key}"] = abs(a[key] / b[key] - 1.0) <= IDENTITY_BAND
        gates[f"sheet {a['tag']} az"] = abs(a["azimuth_deg"] - b["azimuth_deg"]) <= IDENTITY_BAND * 360.0

    built_rod = s10["captured"]["rod"]
    v_space = rod.fields["P1"].e_complex.function_space
    owned, n_wall = _wall_dofs(rod.mesh, built_rod["mesh_facet_tags"], v_space)
    wall = {p: _wall_max_e(rod.fields[p].e_complex, owned, comm) for p in PORT_IDS}
    gates["E"] = max(wall.values()) <= PEC_WALL_MAX_E
    gates["wall_dofs"] = n_wall > 0
    _p0(f"[WF-15 s10] (E) Study's wall dofs {n_wall}, max|E| {wall} (bound {PEC_WALL_MAX_E:g})")

    off = s10["off"]
    wall_rod = _parent_wall_count(rod.mesh, built_rod["mesh_facet_tags"])
    wall_off = _parent_wall_count(off.mesh, s10["captured"]["off"]["mesh_facet_tags"])
    gates["off_cells"] = abs(off.cells / ANS4_EXACT_STRIP_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    gates["wall_rod"] = wall_rod > 0
    gates["wall_off"] = wall_off == 0
    _p0(f"[WF-15 s10] off Study cells {off.cells} (record {ANS4_EXACT_STRIP_CELL_RECORD}, band "
        f"{CELL_COUNT_BAND}); tag-402 owned facets rod {wall_rod}, off {wall_off}")

    for pid in PORT_IDS:  # (PP) printed, never asserted: its first reading on this mesh.
        ex = s10["shares"][pid]
        r = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        _p0(f"[WF-15 s10] (PP) {pid} direct: {r:.3e} (DISCRETE_IDENTITY_RTOL {DISCRETE_IDENTITY_RTOL:g}, "
            f"printed, not asserted)")
    tips = list(TIPS)
    for pid in PORT_IDS:
        row = []
        for i, tip in enumerate(tips):
            a = float(s10["rod_sar"][pid][i]["averaged_sar_w_per_kg"])
            b = float(s10["off_sar"][pid][i]["averaged_sar_w_per_kg"])
            q = a / b if b != 0.0 else float("nan")
            row.append(f"{tip} {q:.6f} (rod {a:.6e}, off {b:.6e} W/kg)")
        _p0(f"[WF-15 s10] tip 1 g Q {pid} (printed, no record on this mesh): " + "; ".join(row))
    red = [k for k, v in gates.items() if not v]
    _p0(f"[WF-15 s10] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red
