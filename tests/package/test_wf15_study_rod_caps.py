"""`WF-15` step 13 (§9 item 802, §10 U121 attempt 2): the rod with capacitor sheets
through the ``Study``, fixture (a) of the 2026-11-03 10:30 review ruling.

The fixture is the PEC rod on U97′'s hole-free exact-rectangle rod build
(``ROD_EXACT_CELL_RECORD`` 168 507, ``tests/package/test_wf15_study_implant_rod.py:67``):
``Study(FourPortBirdcageGeometry(implant_rod=ROD), gel, CapacitorSheets(C_F, esr_ohm=0.2,
on=("P2",)), frequency_hz=64e6, degree=1)`` -- no ``as_hole``, no ``LeontovichWall``, strip
fraction ``GATED_WIDTH_FRACTION`` (the default) on both builds, exact rectangles, κ off.
``C_F`` (``test_ans7_exact_strip``), ``ESR_OHM`` 0.2 and the ``P2`` sheet
(``PLANT_INDEX``) are `PORT-24`'s; its conductor model (copper on the hole mesh) is **not**
used: the conductor substitution copper -> PEC is the review's ruling (1), for the 11-04
weekly to ratify or revert.

Knob ``FEM_EM_WF15_S13`` = ``1``; unset, the coding tests run (no mesh) and the solve node
skips.  Every band is imported; no constant is set here.

Asserted: (A0) the rod Study's cells equal ``ROD_EXACT_CELL_RECORD`` and the rod-free
Study's ``ANS4_EXACT_STRIP_CELL_RECORD``; (Id) on both builds the ``Study``'s 3x3 equals
the direct build (``_build_rod`` + ``build_four_port_sweep(build_only=True)`` +
`PORT-24`'s ``_terminated_three_port``) at ``IDENTITY_BAND``; (i) ``RECIPROCITY_BAND``;
(ii) sigma_max <= 1 + ``PASSIVITY_SIGMA_TOLERANCE``; (E) rod wall <= ``PEC_WALL_MAX_E`` on
the three terminated fields.  Controls (no mesh): the degree-2 rod Study raises its
``ValueError``; the copper-hole rod Study raises ``TypeError``.

Control 2's exception, disclosed: item 802's brief named "the generator's ValueError" (the
item 785 Wcells refusal).  With the lift re-pointed to PEC only (review ruling (3)), the
copper-hole rod Study is refused first by ``_refuse_rod_route``'s own ``TypeError``, before
the generator is reached -- the refusal ``main`` raised before this step, i.e. unchanged.
The target (copper-hole stays refused; lift on PEC only) is kept; the slot owner ruled the
``TypeError`` is what is asserted.

Printed, as a record: ``Study.tune(k=1)``'s ``C_mode`` on the rod PEC four-port beside the
rod-free one (the rod's detuning on PEC).  Leak guard: no S entry is printed, only
relative norms and singular values.  Not claimed: degree 2, copper, a tuned rod-loaded coil.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.ports.circuit import reduce_terminated_ports
from fem_em_solver.study_yaml import StudyConfigError

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD, _build_rod
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf15_study_implant_rod import (
    RECORD_WIDTH,
    ROD_EXACT_CELL_RECORD,
    _load,
    _mapping,
    _phantom,
    _rod_yaml,
)
from tests.validation.test_ans7_exact_strip import C_F
from tests.validation.test_imp3_rod_birdcage import PEC_WALL_MAX_E, _wall_dofs, _wall_max_e
from tests.validation.test_imp3_rod_heating_phantom import _parent_wall_count
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port24_series_rc_sheet import (
    DEGREE,
    ESR_OHM,
    FREQUENCY_HZ,
    PLANT_INDEX,
    _rc,
    _rel,
    _terminated_three_port,
)
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    build_four_port_sweep,
)
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

KNOB = "FEM_EM_WF15_S13"
PORTS = ("P1", "P2", "P3", "P4")
SHEET = PORTS[PLANT_INDEX]

needs_knob = pytest.mark.skipif(os.environ.get(KNOB) != "1", reason=f"{KNOB} unset")
coding = pytest.mark.skipif(os.environ.get(KNOB) == "1", reason=f"{KNOB} set: coding tests skip")
COMM = MPI.COMM_WORLD


def _p0(*args):
    if COMM.rank == 0:
        print(*args, flush=True)


def _rod_geometry():
    return FourPortBirdcageGeometry(implant_rod=ROD)


def _study(geometry, excitation, degree=1, conductor=None):
    return Study(geometry, _phantom(), excitation, frequency_hz=FREQUENCY_HZ, degree=degree,
                 conductor=conductor)


def _caps():
    return CapacitorSheets(C_F, ESR_OHM, on=(SHEET,))


class _Reached(Exception):
    pass


@coding
def test_s13_refusal_lifts_on_the_gated_path_only(monkeypatch):
    """The PEC rod at degree 1 reaches the builder; degree 2 and copper-hole still raise."""
    calls = []

    def _build(self, comm=COMM):
        calls.append(self)
        raise _Reached

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _build)
    assert FREQUENCY_HZ == 64.0e6 and DEGREE == 1 and ESR_OHM == 0.2
    assert _rod_geometry().strip_fraction == GATED_WIDTH_FRACTION
    assert not _rod_geometry().as_hole
    with pytest.raises(_Reached):
        _study(_rod_geometry(), _caps()).solve(COMM)
    # Control 1: the degree-2 rod Study keeps its ValueError.
    with pytest.raises(ValueError, match="degree 2"):
        _study(_rod_geometry(), _caps(), degree=2).solve(COMM)
    # Control 2: the copper-hole rod Study.  The brief named the generator's ValueError;
    # review ruling (3) leaves the Study's own TypeError first (unchanged from `main`), and
    # the slot owner ruled that is what is asserted.  Refused before any mesh.
    copper = FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION,
                                      implant_rod=ROD)
    with pytest.raises(TypeError, match="CapacitorSheets") as err:
        _study(copper, _caps(), conductor=LeontovichWall(study_mod.COPPER_SIGMA)).solve(COMM)
    # The other entry points keep the TypeError on the PEC rod (they never had the route).
    with pytest.raises(TypeError, match="CapacitorSheets"):
        _study(_rod_geometry(), _caps()).solve_sweep((64.0e6, 128.0e6), COMM)
    assert len(calls) == 1
    _p0(f"\n[WF-15 s13] (Ref) PEC rod degree 1 reaches the builder; degree 2 ValueError; "
        f"copper-hole rod TypeError ({err.value}); builds {len(calls)}")


@coding
def test_s13_yaml_admits_the_gated_path_only(tmp_path):
    caps = {"kind": "capacitor_sheets", "c_f": C_F, "esr_ohm": ESR_OHM, "sheets": [SHEET]}
    doc = _mapping(study={"frequency_hz": FREQUENCY_HZ, "degree": 1}, excitation=caps,
                   implant=_rod_yaml())
    st = _load(tmp_path, "gated", doc)
    assert st.geometry == _rod_geometry() and st.conductor is None
    assert isinstance(st.excitation, CapacitorSheets) and st.excitation.on == (SHEET,)
    hole = {"kind": "four_port_birdcage", "preset": "f_small", "as_hole": True}
    refused = {
        "copper hole": ({**doc, "geometry": hole, "conductor": {"kind": "leontovich"}},
                        "implant"),
        "degree 2": ({**doc, "study": {"frequency_hz": FREQUENCY_HZ, "degree": 2}}, "study.degree"),
    }
    for name, (mapping, field) in refused.items():
        with pytest.raises(StudyConfigError) as err:
            _load(tmp_path, name.replace(" ", "_"), mapping)
        _p0(f"[WF-15 s13] YAML {name}: {err.value.field}: {err.value.reason}")
        assert err.value.field == field, (name, err.value.field)


def _direct_three_port(built, label):
    """The direct route on a given PEC build: U97′'s case, `PORT-24`'s P2 3x3."""
    case = build_four_port_sweep(frequency_hz=FREQUENCY_HZ, built=built, build_only=True)
    assert all(s.width_correction_kappa is None for s in case["specs"])
    t0 = time.perf_counter()
    s3 = _terminated_three_port(case["problem"], case["port_defs"], case["specs"],
                                case["facet_tags"], None, PLANT_INDEX, _rc(ESR_OHM))
    _p0(f"[WF-15 s13] direct {label}: 3x3 in {time.perf_counter() - t0:.1f} s")
    return s3, int(case["cells"])


@pytest.fixture(scope="module")
def s13():
    assert COMM.size == RECORD_WIDTH, f"-n {COMM.size}"
    real_build = FourPortBirdcageGeometry.build
    cache = {}
    counts = {"real": 0, "hits": 0}

    def _build(self, comm_=COMM):
        key = "rod" if self.implant_rod else "off"
        if key in cache and cache[key][0] == self:
            counts["hits"] += 1
            return cache[key][1]
        counts["real"] += 1
        out = real_build(self, comm_)
        cache[key] = (self, out)
        return out

    out = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        # The rod-free control first (premise check: never built before this step).
        for key, geo in (("off", FourPortBirdcageGeometry()), ("rod", _rod_geometry())):
            t0 = time.perf_counter()
            res = _study(geo, _caps()).solve(COMM)
            t_caps = time.perf_counter() - t0
            t0 = time.perf_counter()
            st50 = _study(geo, LumpedPortDrive())
            s0 = np.asarray(st50.solve(COMM).s_matrix, dtype=np.complex128)
            tune = st50.tune(k=1, comm=COMM)
            t_50 = time.perf_counter() - t0
            out[key] = {"res": res, "s0": s0, "tune": tune, "built": cache[key][1]}
            _p0(f"\n[WF-15 s13] {key} Study: cells {res.cells}, caps 3x3 {t_caps:.1f} s, "
                f"50 ohm 4x4 + tune {t_50:.1f} s")
    _p0(f"[WF-15 s13] build cache: {counts['real']} real, {counts['hits']} hits")
    assert counts["real"] == 2

    for key, rod in (("off", None), ("rod", ROD)):
        t0 = time.perf_counter()
        built = _build_rod(implant_rod=rod, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
        out[key]["direct"], out[key]["direct_cells"] = _direct_three_port(built, key)
        _p0(f"[WF-15 s13] direct {key}: build {built[4]:.1f} s, total "
            f"{time.perf_counter() - t0:.1f} s")

    # (E): the rod Study's own three terminated fields on its tag-402 wall.
    from fem_em_solver.ports.lumped import run_lumped_sheet_port_case

    res = out["rod"]["res"]
    owned = n_wall = None
    wall = {}
    for pd in res.port_definitions:
        _r, fields = run_lumped_sheet_port_case(
            res.problem, res.all_port_definitions, res.specs, facet_tags=res.facet_tags,
            driven_port_id=pd.port_id, verbose=False, return_fields=True,
            extra_bilinear_terms=res.conductor_terms, degree=res.degree)
        if owned is None:
            owned, n_wall = _wall_dofs(res.mesh, out["rod"]["built"]["mesh_facet_tags"],
                                       fields.e_complex.function_space)
        wall[pd.port_id] = _wall_max_e(fields.e_complex, owned, COMM)
    out["wall"], out["n_wall"] = wall, n_wall
    return out


@complex_only
@needs_knob
def test_s13_all_gates(s13):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    rod, off = s13["rod"], s13["off"]
    gates = {}
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    for key in ("off", "rod"):
        arm = s13[key]
        s_st = np.asarray(arm["res"].s_matrix, dtype=np.complex128)
        rel = _rel(s_st, arm["direct"])
        sig = np.linalg.svd(s_st, compute_uv=False)
        recip = float(np.linalg.norm(s_st - s_st.T) / np.linalg.norm(s_st))
        ident = _rel(s_st, reduce_terminated_ports(arm["s0"], z0, {PLANT_INDEX: _rc(ESR_OHM)}))
        gates[f"{key} Id"] = rel <= IDENTITY_BAND
        gates[f"{key} cells_direct"] = arm["res"].cells == arm["direct_cells"]
        gates[f"{key} i"] = recip <= RECIPROCITY_BAND
        gates[f"{key} ii_sigma"] = float(sig[0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
        gates[f"{key} ports"] = [p.port_id for p in arm["res"].port_definitions] == [
            p for p in PORTS if p != SHEET]
        _p0(f"[WF-15 s13] {key}: cells Study {arm['res'].cells}, direct {arm['direct_cells']}; "
            f"(Id) ||S_study - S_direct||/||S_direct|| {rel:.3e} (IDENTITY_BAND "
            f"{IDENTITY_BAND:.0e}); (i) {recip:.3e} (band {RECIPROCITY_BAND:.0e}); (ii) "
            f"sigma_max {float(sig[0]):.12f} (<= 1 + {PASSIVITY_SIGMA_TOLERANCE:g}); "
            f"reduce(S0) vs Study {ident:.3e} (printed only)")
    gates["A0 rod"] = rod["res"].cells == ROD_EXACT_CELL_RECORD
    gates["A0 off"] = off["res"].cells == ANS4_EXACT_STRIP_CELL_RECORD
    gates["E"] = max(s13["wall"].values()) <= PEC_WALL_MAX_E
    gates["wall_dofs"] = s13["n_wall"] > 0
    wall_rod = _parent_wall_count(rod["res"].mesh, rod["built"]["mesh_facet_tags"])
    wall_off = _parent_wall_count(off["res"].mesh, off["built"]["mesh_facet_tags"])
    gates["wall_rod"] = wall_rod > 0
    gates["wall_off"] = wall_off == 0
    _p0(f"[WF-15 s13] (A0) rod {rod['res'].cells} (record {ROD_EXACT_CELL_RECORD}), off "
        f"{off['res'].cells} (record {ANS4_EXACT_STRIP_CELL_RECORD})")
    _p0(f"[WF-15 s13] (E) wall dofs {s13['n_wall']}, max|E| {s13['wall']} (bound "
        f"{PEC_WALL_MAX_E:g}); tag-402 owned facets rod {wall_rod}, off {wall_off}")
    c_rod, c_off = rod["tune"].c_mode_farad, off["tune"].c_mode_farad
    _p0(f"[WF-15 s13] tune(k=1) C_mode rod {c_rod:.9e} F, rod-free {c_off:.9e} F: "
        f"rod / rod-free - 1 = {c_rod / c_off - 1.0:+.6e} (record, printed only)")
    red = [k for k, v in gates.items() if not v]
    _p0(f"[WF-15 s13] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red
