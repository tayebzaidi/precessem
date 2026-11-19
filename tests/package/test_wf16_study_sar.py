"""`WF-16` step 4 (§9 item 321, chain U23): ``Study.sar(mass_g=10)`` on the exact four-port.

`MAT-4` step 4's C95.3 operator (``post.sar.mass_averaged_sar``) runs through the
``Study`` on item 192's fixture (``FourPortBirdcageGeometry.f_small()``, exact
rectangles, κ off) at 64 MHz, four single drives.  Asserted: (R) the
field-keeping solve reproduces step 1's S record at ``LEG_D0_REPRODUCTION_BAND``;
(i) the whole-phantom coverage identity; (ii) the 10 g C4 quadrant identity at
``C4_COVARIANCE_BAND``; (iii) ``Study.sar`` == the direct operator call; (iv)
containment; (c1) the mis-paired far-side ball misses the band (sign asserted,
size printed); (c2) refusals.  Every band is imported; none is set here.

Printed, never asserted: the 4 × 1 diagonal and the 10 g maximum over the four
evaluated centres at 1 V drive (NOT a spatial search, NOT a C95.3 compliance
figure), the kernel mass, the phantom cell count and ``a_10g / h``.

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
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
    mean_sar,
)
from fem_em_solver.study import (
    SAR_QUADRATURE_DEGREE,
    QuasiStaticSphereField,
    SphereInBoxGeometry,
)
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_phantom_resolution import _global_tag_cell_count
from tests.package.test_wf16_study_four_port import (
    DEGREE2_RECORD_WIDTH,
    _degree,
    _p0,
    _record,
    _rel,
    needs_degree,
)
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND
from tests.validation.test_birdcage_b1_quadrature import (
    QUADRATURE_STEP_DEG,
    _port_index,
)
from tests.validation.test_birdcage_sar_mass_averaged import (
    CENTRE_RADIUS_M,
    CONTROL_CEILING,
    CONTROL_PREDICTED,
    EXACT_IDENTITY_RTOL,
    WHOLE_PHANTOM_BALL_RADIUS_M,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_mass_averaged_sar_standard_masses import (
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

SAR = "averaged_sar_w_per_kg"


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


def test_wf16_sar_refusals_and_pins():
    """(c2) refusals before any mesh; the pinned quadrature degree is MAT-4's."""
    assert SAR_QUADRATURE_DEGREE == QUADRATURE_DEGREE
    sphere = Study(
        SphereInBoxGeometry(0.05, 0.15, 0.01, 0.03),
        _phantom(),
        QuasiStaticSphereField(),
        frequency_hz=FREQUENCY_64_HZ,
    )
    with pytest.raises(TypeError):
        sphere.sar(centres=[(0.0, 0.0, 0.0)])
    study = _port_study()
    with pytest.raises(ValueError):
        study.sar(mass_g=0, centres=[(0.0, 0.0, 0.0)])
    with pytest.raises(ValueError):
        study.sar(mass_g=10.0, radius_m=0.01, centres=[(0.0, 0.0, 0.0)])
    with pytest.raises(ValueError):
        study.sar(mass_g=None, centres=[(0.0, 0.0, 0.0)])
    assert getattr(study, "fields_result", None) is None, "a refusal built a mesh"


def test_wf16_sar_containment():
    """(iv) on the Study's own geometry fields (no solve)."""
    g = FourPortBirdcageGeometry.f_small()
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=_port_study().density_kg_m3)
    assert CENTRE_RADIUS_M + a10 < g.phantom_radius
    circum = float(np.hypot(g.phantom_radius, 0.5 * g.phantom_height))
    nearest_conductor = min(g.ring_radius - 0.5 * g.leg_width, g.ring_radius - g.ring_minor_radius)
    assert circum < WHOLE_PHANTOM_BALL_RADIUS_M < nearest_conductor
    _p0(
        f"\n[WF-16 s4] (iv) r0 + a_10g = {(CENTRE_RADIUS_M + a10) * 1e3:.4f} mm < "
        f"{g.phantom_radius * 1e3:.1f} mm; circumradius {circum:.4f} m < "
        f"{WHOLE_PHANTOM_BALL_RADIUS_M} m < conductor {nearest_conductor:.4f} m"
    )


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
    assert study.result is None, "the field-keeping solve wrote self.result"
    assert set(res.fields) == {"P1", "P2", "P3", "P4"}

    sheets = {f"P{s['tag'] - study.geometry.sheet_iface}": s for s in res.sheets}
    az = {pid: float(s["azimuth_deg"]) for pid, s in sheets.items()}
    idx = {pid: _port_index(az[pid], az["P1"], QUADRATURE_STEP_DEG) for pid in az}
    by_k = {k: pid for pid, k in idx.items()}
    assert sorted(by_k) == [0, 1, 2, 3]
    centres = {
        k: (
            CENTRE_RADIUS_M * float(np.cos(np.radians(az[by_k[k]]))),
            CENTRE_RADIUS_M * float(np.sin(np.radians(az[by_k[k]]))),
            0.0,
        )
        for k in range(4)
    }

    t0 = time.perf_counter()
    diag = {
        k: study.sar(mass_g=10.0, centres=[centres[k]], drives=[by_k[k]], comm=comm)[by_k[k]][0]
        for k in range(4)
    }
    control = {
        k: study.sar(mass_g=10.0, centres=[centres[(k + 2) % 4]], drives=[by_k[k]], comm=comm)[
            by_k[k]
        ][0]
        for k in range(4)
    }
    whole = study.sar(
        radius_m=WHOLE_PHANTOM_BALL_RADIUS_M,
        centres=[(0.0, 0.0, 0.0)],
        drives=["P1"],
        comm=comm,
    )["P1"][0]
    t_ops = time.perf_counter() - t0

    # (iii) and (i)'s comparators, built here from the test's own ρ and radius.
    f1 = res.fields["P1"]
    all_tags = global_cell_tag_set(res.mesh, res.cell_tags)
    rho = build_density_field(
        res.mesh,
        study.density_kg_m3,
        cell_tags=res.cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != study.geometry.phantom_tag},
    )
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=study.density_kg_m3)
    direct = mass_averaged_sar(
        f1.e_complex,
        sigma=f1.sigma_field,
        rho=rho,
        center=centres[idx["P1"]],
        radius=a10,
        comm=comm,
        quadrature_degree=QUADRATURE_DEGREE,
    )
    tagged = mean_sar(
        f1.e_complex,
        sigma=f1.sigma_field,
        rho=study.density_kg_m3,
        cell_tags=res.cell_tags,
        subdomain_ids=study.geometry.phantom_tag,
        comm=comm,
    )
    phantom_cells = _global_tag_cell_count(
        res.mesh, res.cell_tags, study.geometry.phantom_tag, comm
    )
    _p0(
        f"\n[WF-16 s4] degree {degree} -n {comm.size}: field-keeping solve "
        f"{t_solve:.2f} s (mesh {res.mesh_time:.2f} s, sweep {res.sweep_time:.2f} s); "
        f"9 operator calls {t_ops:.2f} s at quadrature_degree {SAR_QUADRATURE_DEGREE}; "
        f"{res.cells} cells, {phantom_cells} phantom (tag {study.geometry.phantom_tag}) "
        f"cells; resolution {study.geometry.resolution} m, a_10g = {a10 * 1e3:.4f} mm, "
        f"a_10g / h = {a10 / study.geometry.resolution:.4f}"
    )
    report_peak_rss(comm, label="[WF-16 s4] summed peak ru_maxrss")
    return {
        "degree": degree,
        "study": study,
        "res": res,
        "by_k": by_k,
        "idx": idx,
        "centres": centres,
        "diag": diag,
        "control": control,
        "whole": whole,
        "direct": direct,
        "tagged": tagged,
        "a10": a10,
    }


@complex_only
@needs_degree
def test_wf16_sar_r_same_fixture(case):
    """(R) keeping fields did not move the solve: S vs step 1's record."""
    rel = _rel(case["res"].s_matrix, _record(case["degree"]))
    _p0(
        f"\n[WF-16 s4] (R) ||S_keep - record|| / ||record|| = {rel:.3e} "
        f"(ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    assert rel <= LEG_D0_REPRODUCTION_BAND


@complex_only
@needs_degree
def test_wf16_sar_i_whole_phantom_coverage(case):
    """(i) whole ball == tagged mean_sar power; mass == ρ·V_phantom."""
    w, t = case["whole"], case["tagged"]
    rho = case["study"].density_kg_m3
    d_p = abs(w["dissipated_power_w"] / t["dissipated_power_w"] - 1.0)
    d_m = abs(w["mass_kg"] / (rho * t["volume_m3"]) - 1.0)
    _p0(
        f"\n[WF-16 s4] (i) ball power {w['dissipated_power_w']:.12e} W, tagged "
        f"{t['dissipated_power_w']:.12e} W: rel {d_p:.3e}; mass {w['mass_kg']:.12e} kg "
        f"vs rho*V {rho * t['volume_m3']:.12e} kg: rel {d_m:.3e} "
        f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e})"
    )
    assert d_p <= EXACT_IDENTITY_RTOL
    assert d_m <= EXACT_IDENTITY_RTOL


@complex_only
@needs_degree
def test_wf16_sar_ii_ten_gram_c4(case):
    """(ii) the 10 g C4 quadrant identity, the gate."""
    diag, by_k = case["diag"], case["by_k"]
    pairs = {k: abs(diag[(k + 1) % 4][SAR] - diag[k][SAR]) / diag[k][SAR] for k in range(4)}
    lines = "\n".join(
        f"    k={k} ({by_k[k]} at c_{k}): SAR_10g = {diag[k][SAR]:.9e} W/kg, mass "
        f"{diag[k]['mass_kg']:.9e} kg vs rho*4/3*pi*a^3 "
        f"{case['study'].density_kg_m3 * 4.0 / 3.0 * np.pi * case['a10'] ** 3:.9e} kg"
        for k in range(4)
    )
    k_max = max(range(4), key=lambda k: diag[k][SAR])
    _p0(
        f"\n[WF-16 s4] (ii) 10 g diagonal (printed):\n{lines}\n"
        "    cyclic pairs |SAR(c_k+1;k+1) - SAR(c_k;k)| / SAR(c_k;k): "
        + ", ".join(f"{k}->{(k + 1) % 4} {pairs[k] * 100:.4f}%" for k in range(4))
        + f" (ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%)\n"
        f"    10 g max over 4 evaluated centres at 1 V drive, NOT a spatial search, "
        f"NOT a C95.3 compliance figure: {diag[k_max][SAR]:.9e} W/kg ({by_k[k_max]})"
    )
    for k in range(4):
        assert pairs[k] <= C4_COVARIANCE_BAND, (k, pairs[k])


@complex_only
@needs_degree
def test_wf16_sar_iii_api_equals_operator(case):
    """(iii) Study.sar dict == the direct operator call on the kept field."""
    api = case["diag"][case["idx"]["P1"]]
    direct = case["direct"]
    worst = 0.0
    for key, value in direct.items():
        if isinstance(value, (int, float)) and value != 0.0:
            worst = max(worst, abs(api[key] - value) / abs(value))
    _p0(
        f"\n[WF-16 s4] (iii) Study.sar vs mass_averaged_sar: max rel {worst:.3e} "
        f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e}); keys {sorted(direct)}"
    )
    assert set(api) == set(direct)
    assert worst <= EXACT_IDENTITY_RTOL


@complex_only
@needs_degree
def test_wf16_sar_c1_mispaired_control(case):
    """(c1) far-side ball: sign asserted > band, size printed (predicted only)."""
    diag, control = case["diag"], case["control"]
    miss = {k: abs(control[k][SAR] - diag[k][SAR]) / diag[k][SAR] for k in range(4)}
    _p0(
        "\n[WF-16 s4] (c1) mis-paired far-side ball |SAR(c_k+2;k) - SAR(c_k;k)| / "
        "SAR(c_k;k): "
        + ", ".join(f"k={k} {miss[k] * 100:.4f}%" for k in range(4))
        + f" (sign ASSERTED > {C4_COVARIANCE_BAND * 100:.1f}%; size PREDICTED "
        f"≈ {CONTROL_PREDICTED * 100:.0f}%, ceiling {CONTROL_CEILING * 100:.0f}%, "
        "75-95% at 64 MHz, not asserted)"
    )
    for k in range(4):
        assert miss[k] > C4_COVARIANCE_BAND, (k, miss[k])
