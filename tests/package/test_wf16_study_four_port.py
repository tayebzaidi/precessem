"""`WF-16` step 1 (§9 item 244, chain U4): the ``Study`` on the F-small 4-port birdcage.

The same ``Study`` class `WF-15` landed builds item 192's ×1 fixture on exact
port rectangles with κ off (``FourPortBirdcageGeometry.f_small()``), drives it
with ``LumpedPortDrive``, solves at 64 MHz and returns S and a Touchstone
file.  Both are held to `ANS-4` step 5's registered exact-rectangle records at
``LEG_D0_REPRODUCTION_BAND`` (1e-9), at the record's own degree and width:
degree 1 at ``-n 2`` (``ANS4_EXACT_STRIP_S_P1_RECORDS``), degree 2 — the Study's
default — at ``-n 8`` (``ANS4_EXACT_STRIP_S_P2_X1``).  `PORT-11`'s gates are
imported.  Every band is imported; none is set here.

Knob ``FEM_EM_WF16_DEGREE`` in {1, 2}; unset, the solve nodes skip; any other
value raises.  The width guard is asserted, never skipped.
"""

from __future__ import annotations

import dataclasses
import os
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.core import build_material_fields, build_mu_r_field
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.utils.instrumentation import report_peak_rss

from tests.complex_mode import complex_only
from tests.materials.test_phantom_material_model import _dg0_owned_diff
from tests.mesh.test_birdcage_leg_gaps import LEG_GAP_LENGTH
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND, SHEET_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import (
    CELL_COUNT_BAND,
    CONDUCTOR_RESOLUTION,
)
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import (
    AIR_PADDING,
    COIL_LENGTH,
    LEG_COUNT,
    LEG_SPACING,
    LEG_WIDTH,
    PHANTOM_HEIGHT,
    PHANTOM_RADIUS,
    PORT_BOX_SIZE,
    RESOLUTION,
    RING_MINOR_RADIUS,
    RING_RADIUS,
)
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    ANS4_EXACT_STRIP_S_P1_RECORDS,
    ANS4_EXACT_STRIP_S_P2_X1,
    SPEC_AREA_M2,
    SPEC_HEIGHT_M,
    SPEC_WIDTH_M,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import (
    ANS4_EXACT_STRIP_CELL_RECORD,
    STRIP_BAND,
)
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    TERMINATED_PORT_IMPEDANCE_OHM,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_leg_offset_sweep import _four_port_rung
from tests.validation.test_port_birdcage_lumped_column import (
    ADJACENT_SPREAD_BAND,
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "logs" / "wf16"
DEGREE_ENV = "FEM_EM_WF16_DEGREE"
DEGREE2_RECORD_WIDTH = 8  # `P2_X1`'s version tag (`-n 8`), not a band.


def _degree():
    raw = os.environ.get(DEGREE_ENV)
    if raw is None or raw == "":
        return None
    if raw not in ("1", "2"):
        raise ValueError(f"{DEGREE_ENV}={raw!r}: only 1 or 2")
    return int(raw)


needs_degree = pytest.mark.skipif(_degree() is None, reason=f"{DEGREE_ENV} unset")


def _record(degree):
    table = ANS4_EXACT_STRIP_S_P1_RECORDS if degree == 1 else ANS4_EXACT_STRIP_S_P2_X1
    return np.asarray(table["64"], dtype=np.complex128)


def _rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _fields_of(obj):
    if dataclasses.is_dataclass(obj):
        return {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj)}
    return dict(vars(obj))


@pytest.fixture(scope="module")
def leg():
    degree = _degree()
    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH if degree == 1 else DEGREE2_RECORD_WIDTH
    # Width guard, asserted (a record holds at its own width, `OPS-41`).
    assert comm.size == width, (
        f"degree {degree} record is asserted at -n {width}; got -n {comm.size}"
    )
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    kwargs = {"frequency_hz": FREQUENCY_64_HZ}
    if degree == 1:
        kwargs["degree"] = 1  # degree 2 is the Study's default: no argument.
    study = Study(
        FourPortBirdcageGeometry.f_small(), phantom, LumpedPortDrive(), **kwargs
    )
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_study = time.perf_counter() - t0
    path = OUT_DIR / f"wf16_f_small_p{degree}_n{comm.size}.s4p"
    t0 = time.perf_counter()
    written = res.write_touchstone(path)
    t_ts = time.perf_counter() - t0
    t0 = time.perf_counter()
    direct = _four_port_rung(
        "wf16-direct",
        np.zeros(4),
        frequency_hz=FREQUENCY_64_HZ,
        degree=degree,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    t_direct = time.perf_counter() - t0
    _p0(
        f"\n[WF-16] degree {degree} -n {comm.size}: Study {t_study:.2f} s (mesh "
        f"{res.mesh_time:.2f} s, sweep {res.sweep_time:.2f} s); Touchstone "
        f"{t_ts:.3f} s; direct call {t_direct:.2f} s (sweep "
        f"{direct['sweep_time']:.2f} s)"
    )
    report_peak_rss(comm, label="[WF-16] summed peak ru_maxrss")
    return {
        "degree": degree,
        "study": study,
        "res": res,
        "path": written,
        "direct": direct,
        "phantom": phantom,
    }


def test_wf16_geometry_constants_and_refusals():
    """(K) constants restated == imported; κ and ragged refused; default degree 2."""
    g = FourPortBirdcageGeometry.f_small()
    expected = {
        "leg_count": LEG_COUNT,
        "ring_radius": RING_RADIUS,
        "leg_width": LEG_WIDTH,
        "leg_spacing": LEG_SPACING,
        "coil_length": COIL_LENGTH,
        "ring_minor_radius": RING_MINOR_RADIUS,
        "phantom_radius": PHANTOM_RADIUS,
        "phantom_height": PHANTOM_HEIGHT,
        "port_box_size": PORT_BOX_SIZE,
        "leg_gap_length": LEG_GAP_LENGTH,
        "air_padding": AIR_PADDING,
        "resolution": RESOLUTION,
        "conductor_resolution": CONDUCTOR_RESOLUTION,
        "strip_fraction": GATED_WIDTH_FRACTION,
        "sheet_iface": SHEET_IFACE,
        "port_lower": PORT_LOWER,
        "port_upper": PORT_UPPER,
        "conductor_tag": CONDUCTOR_CELL_TAG,
        "phantom_tag": PHANTOM_CELL_TAG,
        "conductor_sigma_s_per_m": SIGMA_WIRE_S_PER_M,
        # `WF-16` step 7 (§9 item 324): the additive hole-route field; the
        # default stays the solid conductor (item 192's fixture).
        "as_hole": False,
        # `WF-15` step 10 (§9 item 709): additive; the default stays no implant.
        "implant_rod": None,
    }
    for name, value in expected.items():
        assert getattr(g, name) == value, name
    assert set(_fields_of(g)) == set(expected)
    with pytest.raises(TypeError):
        FourPortBirdcageGeometry(width_correction_kappa=1.0)
    with pytest.raises(ValueError):
        FourPortBirdcageGeometry(strip_fraction=None)
    drive = LumpedPortDrive()
    assert drive.z0_ohm == REFERENCE_IMPEDANCE_OHM
    assert drive.port_impedance_ohm == TERMINATED_PORT_IMPEDANCE_OHM
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    study = Study(g, phantom, drive, frequency_hz=FREQUENCY_64_HZ)
    assert study.degree == 2
    with pytest.raises(TypeError):
        study.mean_sar()


@complex_only
@needs_degree
def test_wf16_upstream_identity(leg):
    """(F) SPEC rectangle by CAD identity; (U) Study builds the direct problem; (c3)."""
    comm = MPI.COMM_WORLD
    res, direct = leg["res"], leg["direct"]
    # (F)
    ratio = res.cells / ANS4_EXACT_STRIP_CELL_RECORD
    _p0(f"[WF-16] (F) cells {res.cells} / record {ANS4_EXACT_STRIP_CELL_RECORD}")
    assert abs(ratio - 1.0) <= CELL_COUNT_BAND
    for s in res.sheets:
        _p0(
            f"    sheet {s['tag']}: {s['facets']} facets area {s['area']:.12e} "
            f"w {s['w']:.12e} h {s['h']:.12e} centre {s['centre_radial']:.12e}"
        )
        assert abs(s["area"] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["w"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["w_bbox"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["h"] / SPEC_HEIGHT_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["centre_radial"] / RING_RADIUS - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["area_ratio_full"] - 1.0) < SHEET_AREA_BAND, s["tag"]
    # (U) cells, per-port geometry, specs, port definitions.
    assert res.cells == direct["cells"]
    for a, b in zip(res.sheets, direct["sheets"]):
        assert a["tag"] == b["tag"] and a["facets"] == b["facets"]
        for key in ("h", "w", "area"):
            assert abs(a[key] / b[key] - 1.0) <= IDENTITY_BAND, (a["tag"], key)
        assert abs(a["azimuth_deg"] - b["azimuth_deg"]) <= IDENTITY_BAND * 360.0
    assert len(res.specs) == len(direct["specs"]) == 4
    for a, b in zip(res.specs, direct["specs"]):
        fa, fb = _fields_of(a), _fields_of(b)
        assert set(fa) == set(fb)
        assert fa.get("width_correction_kappa") is None
        for key in fa:
            if key in ("gap_height_m", "sheet_width_m"):
                assert abs(fa[key] / fb[key] - 1.0) <= IDENTITY_BAND, key
            else:
                assert fa[key] == fb[key], key
    for a, b in zip(res.port_definitions, direct["port_defs"]):
        assert _fields_of(a) == _fields_of(b)
    assert res.sweep.is_placeholder is False
    # (U) DG0 fields both ways on the Study's mesh.
    p_dir = direct["problem"]

    def _both(phantom):
        ps = res.problem
        s_sig, s_eps = build_material_fields(
            res.mesh,
            ps.material,
            cell_tags=res.cell_tags,
            material_map=ps.material_map,
            phantom_material=phantom,
            phantom_tag=ps.phantom_tag,
        )
        s_mu = build_mu_r_field(
            res.mesh, ps.material, cell_tags=res.cell_tags, material_map=ps.material_map
        )
        d_sig, d_eps = build_material_fields(
            res.mesh,
            p_dir.material,
            cell_tags=res.cell_tags,
            material_map=p_dir.material_map,
        )
        d_mu = build_mu_r_field(
            res.mesh, p_dir.material, cell_tags=res.cell_tags, material_map=p_dir.material_map
        )
        return (
            _dg0_owned_diff(s_sig, d_sig, comm),
            _dg0_owned_diff(s_eps, d_eps, comm),
            _dg0_owned_diff(s_mu, d_mu, comm),
        )

    sig, eps, mu = _both(res.problem.phantom_material)
    _p0(f"[WF-16] (U) DG0 max|d| sigma {sig}, eps_r {eps}, mu_r {mu}")
    assert sig[0] == 0.0 and eps[0] == 0.0 and mu[0] == 0.0
    # (c3) 1-ULP sigma through the phantom slot: exactly 2^-53 on > 0 dofs.
    ulp = GelledSalinePhantomMaterial(np.nextafter(0.5, 1), SALINE_EPSILON_R, FREQUENCY_64_HZ)
    sig_c, _eps_c, _mu_c = _both(ulp)
    _p0(f"[WF-16] (c3) 1-ULP sigma: max|d| {sig_c[0]!r} on {sig_c[1]} owned dofs")
    assert sig_c[0] == 2.0**-53 and sig_c[1] > 0


@complex_only
@needs_degree
def test_wf16_s_equals_record(leg):
    """(R) S equals the record at its degree and width; (c1) cross-degree swap."""
    d = leg["degree"]
    s = leg["res"].s_matrix
    rel = _rel(s, _record(d))
    other = _rel(s, _record(3 - d))
    records = _rel(_record(1), _record(2))
    _p0(
        f"[WF-16] (R{d}) ||S - S_rec||/||S_rec|| {rel:.3e} (band "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}); (c1) vs the other degree's record "
        f"{other:.3e}; record-to-record {records:.4f}"
    )
    assert rel <= LEG_D0_REPRODUCTION_BAND
    assert other > LEG_D0_REPRODUCTION_BAND
    # The brief's hand arithmetic: ≈ 0.049 to two digits, else mis-wired (stop).
    assert 0.0485 <= records < 0.0495, "comparator mis-wired"
    # Printed: against full-precision JSON if present.
    js = REPO / "logs" / "ops61" / f"ours_f05x_p{d}.json"
    _p0(f"[WF-16] full-precision JSON {js.name}: {'present' if js.exists() else 'absent'}")


@complex_only
@needs_degree
def test_wf16_touchstone_equals_record(leg):
    """(T1) file == in-memory S (both readers); (T2) file vs record; (c2) P1<->P2."""
    import skrf

    d = leg["degree"]
    s = leg["res"].s_matrix
    freqs, mats, ports = load_touchstone(leg["path"])
    loaded = np.asarray(mats[0], dtype=np.complex128)
    ntw = skrf.Network(str(leg["path"]))
    sk = np.asarray(ntw.s[0], dtype=np.complex128)
    d_load = float(np.max(np.abs(loaded - s)))
    d_skrf = float(np.max(np.abs(sk - s)))
    _p0(f"[WF-16] (T1) max|dS| load_touchstone {d_load:.3e}, skrf {d_skrf:.3e}")
    assert d_load <= ANCHOR_TOL and d_skrf <= ANCHOR_TOL
    assert abs(float(freqs[0]) - FREQUENCY_64_HZ) <= FREQ_ATOL_HZ
    assert abs(float(ntw.f[0]) - FREQUENCY_64_HZ) <= FREQ_ATOL_HZ
    assert ntw.nports == 4 and tuple(ports) == ("P1", "P2", "P3", "P4")
    assert np.all(np.asarray(ntw.z0) == 50.0)
    rec = _record(d)
    t2 = _rel(loaded, rec)
    perm = [1, 0, 2, 3]
    swapped = loaded[np.ix_(perm, perm)]
    c2 = _rel(swapped, rec)
    rot = [1, 2, 3, 0]
    c_rot = _rel(loaded[np.ix_(rot, rot)], rec)
    c_tr = _rel(loaded.T, rec)
    _p0(
        f"[WF-16] (T2) ||S_file - S_rec||/||S_rec|| {t2:.3e}; (c2) P1<->P2 swap "
        f"{c2:.4f}; printed: cyclic rotation {c_rot:.3e}, transpose {c_tr:.3e}"
    )
    assert t2 <= LEG_D0_REPRODUCTION_BAND
    assert c2 > LEG_D0_REPRODUCTION_BAND


@complex_only
@needs_degree
def test_wf16_port11_gates(leg):
    """(G) `PORT-11` (i), (ii), (iii′) on the Study's own S and Z."""
    s = leg["res"].s_matrix
    z = leg["res"].z_matrix
    recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
    sigma = np.linalg.svd(s, compute_uv=False)
    power = np.sum(np.abs(s) ** 2, axis=0)
    classes = _circulant_classes(z)
    spreads = {n: _class_spread(v) for n, v in classes.items()}
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    worst = max(spreads.values())
    separation = pooled / worst if worst > 0.0 else np.inf
    _p0(
        f"[WF-16] (G) reciprocity {recip:.3e}; sigma_max {float(np.max(sigma)):.12f}; "
        f"max column power {float(np.max(power)):.12f}; spreads "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
        + f"; separation {separation:.2f}x"
    )
    assert recip <= RECIPROCITY_BAND
    assert float(np.max(sigma)) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for v in power:
        assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for name, v in spreads.items():
        assert v <= ADJACENT_SPREAD_BAND, name
    assert separation >= POOLED_SEPARATION_FLOOR


@complex_only
@needs_degree
def test_wf16_study_equals_direct_call(leg):
    """(I) Study S vs the direct call's S, same process."""
    rel = _rel(leg["res"].s_matrix, leg["direct"]["s"])
    _p0(f"[WF-16] (I) ||S_study - S_direct||/||S_direct|| {rel:.3e}")
    assert rel <= LEG_D0_REPRODUCTION_BAND
    assert leg["res"].degree == leg["degree"]


@complex_only
@needs_degree
def test_wf16_record_comparators_are_not_blind(leg):
    """(c1)/(c2) restated as one node: each comparator exceeds the band."""
    d = leg["degree"]
    s = leg["res"].s_matrix
    rec = _record(d)
    c1 = _rel(s, _record(3 - d))
    perm = [1, 0, 2, 3]
    c2 = _rel(s[np.ix_(perm, perm)], rec)
    _p0(f"[WF-16] (c1) {c1:.4f}  (c2) {c2:.4f}  (band {LEG_D0_REPRODUCTION_BAND:.0e})")
    assert c1 > LEG_D0_REPRODUCTION_BAND
    assert c2 > LEG_D0_REPRODUCTION_BAND
