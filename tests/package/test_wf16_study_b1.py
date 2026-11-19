"""`WF-16` step 5 (§9 item 322, chain U24): ``Study.b1_plus()`` on the exact four-port.

`WF-6` step 1d / 2b's gated ``|B₁⁺|`` C4 covariance identity (base plus three
rotated copies, CG1 estimator) runs through the ``Study`` on item 192's fixture
(``FourPortBirdcageGeometry.f_small()``, exact rectangles, κ off) at 64 MHz.
Asserted: (B1) P2@+Δ, P4@−Δ, P3@2Δ against P1 within ``C4_COVARIANCE_BAND``;
(B2) ``Study.b1_plus_at`` == the hand-built CG1 chain on the same kept field;
(B3) every sample point valid on all four drives; the mis-rotation control
(P2 at −Δ) outside the band (sign asserted, size printed).  Every band and
the sample set are imported; none is set here.

Printed, never asserted: ``|B₁⁺|`` at the phantom centre per drive at 1 V (a
record, not homogeneity, not absolute).

Knob ``FEM_EM_WF16_DEGREE`` (step 1's) in {1, 2}; unset, the solve nodes skip.
Width guard asserted: degree 1 at ``-n 2``, degree 2 at ``-n 8``.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.post.faraday import magnetic_flux_density_from_e, project_to_cg1
from fem_em_solver.study import QuasiStaticSphereField, SphereInBoxGeometry
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.package.test_wf16_study_four_port import (
    DEGREE2_RECORD_WIDTH,
    _degree,
    _p0,
    needs_degree,
)
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    PHANTOM_CELL_TAG,
    _read_b1_plus_cg1,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)


def _port_study(degree=2):
    return Study(
        FourPortBirdcageGeometry.f_small(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=degree,
    )


def test_wf16_b1_refusal_and_tag():
    """Sphere Study refused before any mesh; the sample set's tag is the Study's."""
    sphere = Study(
        SphereInBoxGeometry(0.05, 0.15, 0.01, 0.03),
        _phantom(),
        QuasiStaticSphereField(),
        frequency_hz=FREQUENCY_64_HZ,
    )
    with pytest.raises(TypeError):
        sphere.b1_plus()
    with pytest.raises(TypeError):
        sphere.b1_plus_at(np.zeros((1, 3)))
    assert PHANTOM_CELL_TAG == FourPortBirdcageGeometry.f_small().phantom_tag == 3


@pytest.fixture(scope="module")
def case():
    degree = _degree()
    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH if degree == 1 else DEGREE2_RECORD_WIDTH
    assert (
        comm.size == width
    ), f"degree {degree} record is asserted at -n {width}; got -n {comm.size}"
    study = _port_study(degree)
    t0 = time.perf_counter()
    res = study._solve_port_fields(comm)
    t_solve = time.perf_counter() - t0
    assert set(res.fields) == {"P1", "P2", "P3", "P4"}

    az = {f"P{s['tag'] - study.geometry.sheet_iface}": float(s["azimuth_deg"]) for s in res.sheets}
    # Read off the fixture's sheets, never chosen.
    delta_deg = (az["P2"] - az["P1"]) % 360.0
    delta = np.radians(delta_deg)
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})

    t0 = time.perf_counter()
    maps = study.b1_plus(comm=comm)
    t_proj = time.perf_counter() - t0
    images = {
        "P1@0": ("P1", points),
        "P2@+D": ("P2", _rotate_z(points, delta)),
        "P4@-D": ("P4", _rotate_z(points, -delta)),
        "P3@2D": ("P3", _rotate_z(points, 2.0 * delta)),
        "P2@-D (control)": ("P2", _rotate_z(points, -delta)),
    }
    t0 = time.perf_counter()
    read = {
        label: study.b1_plus_at(pts, drives=[pid], comm=comm)[pid]
        for label, (pid, pts) in images.items()
    }
    centre = {
        pid: float(study.b1_plus_at(np.zeros((1, 3)), drives=[pid], comm=comm)[pid][0][0])
        for pid in ("P1", "P2", "P3", "P4")
    }
    t_eval = time.perf_counter() - t0

    # (B2)'s comparator: the hand-built chain on the same kept field.
    hand = project_to_cg1(magnetic_flux_density_from_e(res.fields["P1"].e_complex, study.omega))
    hand_vals, hand_valid = _read_b1_plus_cg1(hand, points)

    _p0(
        f"\n[WF-16 s5] degree {degree} -n {comm.size}: field-keeping solve "
        f"{t_solve:.2f} s (mesh {res.mesh_time:.2f} s, sweep {res.sweep_time:.2f} s); "
        f"4 CG1 projections {t_proj:.2f} s; evaluations {t_eval:.2f} s; "
        f"{res.cells} cells; drive rotation {delta_deg:.6f} deg; "
        f"{points.shape[0]} sample points"
    )
    report_peak_rss(comm, label="[WF-16 s5] summed peak ru_maxrss")
    return {
        "degree": degree,
        "study": study,
        "maps": maps,
        "points": points,
        "read": read,
        "centre": centre,
        "hand": (hand_vals, hand_valid),
    }


@complex_only
@needs_degree
def test_wf16_b1_b3_valid_point_floor(case):
    """(B3) every sample point valid on all four drives and every image."""
    n = case["points"].shape[0]
    counts = {label: int(v.sum()) for label, (_, v) in case["read"].items()}
    _p0(f"\n[WF-16 s5] (B3) valid points per image (of {n}): {counts} (ASSERTED full)")
    assert n > 0
    for label, c in counts.items():
        assert c == n, (label, c, n)


@complex_only
@needs_degree
def test_wf16_b1_b1_c4_covariance(case):
    """(B1) the CG1 C4 covariance identity, the gate."""
    read = case["read"]
    mask = np.logical_and.reduce([v for _, v in read.values()])
    ref = read["P1@0"][0]
    table = {
        label: _relative_l2(vals, ref, mask)
        for label, (vals, _) in read.items()
        if label not in ("P1@0", "P2@-D (control)")
    }
    _p0(
        "\n[WF-16 s5] (B1) ||B1+(image) - B1+(P1)||_2 / ||B1+(P1)||_2: "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in table.items())
        + f" (ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%)\n"
        "    |B1+| at the phantom centre per 1 V drive (printed, not absolute): "
        + ", ".join(f"{p} {v:.9e} T" for p, v in case["centre"].items())
    )
    for label, value in table.items():
        assert value <= C4_COVARIANCE_BAND, (label, value)


@complex_only
@needs_degree
def test_wf16_b1_b2_api_equals_package(case):
    """(B2) Study.b1_plus_at P1 == hand-built magnetic_flux_density_from_e → project_to_cg1."""
    api, api_valid = case["read"]["P1@0"]
    hand, hand_valid = case["hand"]
    assert np.array_equal(api_valid, hand_valid)
    rel = float(
        np.max(np.abs(api[api_valid] - hand[hand_valid])) / np.max(np.abs(hand[hand_valid]))
    )
    _p0(
        f"\n[WF-16 s5] (B2) max|API - hand| / max|B1+| = {rel:.3e} "
        f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e})"
    )
    assert rel <= EXACT_IDENTITY_RTOL


@complex_only
@needs_degree
def test_wf16_b1_control_mis_rotation(case):
    """Control: P2 at the counter-rotated points misses the band (sign asserted)."""
    read = case["read"]
    mask = np.logical_and.reduce([v for _, v in read.values()])
    miss = _relative_l2(read["P2@-D (control)"][0], read["P1@0"][0], mask)
    _p0(
        f"\n[WF-16 s5] control P2@-D vs P1: {miss * 100:.4f}% (sign ASSERTED > "
        f"{C4_COVARIANCE_BAND * 100:.1f}%; size PREDICTED 80-100%, not asserted)"
    )
    assert miss > C4_COVARIANCE_BAND
