"""`WF-16` step 7 (§9 item 324, chain U26): the port ``Study`` on copper.

``FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)``
plus ``Study(..., conductor=LeontovichWall(COPPER_SIGMA))`` reproduces `ANS-6`
step 4's exact-rectangle hole route (§9 item 193) at degree 1, 64 MHz, ``-n 2``:

- (R) the Cu class means equal ``ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]``
  at ``LEG_D0_REPRODUCTION_BAND``;
- (b) ``p_surf > 0`` and the discrete power-identity residual
  <= ``DISCRETE_IDENTITY_RTOL``;
- (g) `TH-14`'s gates (reciprocity, passivity, adjacent spread), imported.

Negative control (rule (e)): the same Study with ``conductor=None`` (the hole
wall left natural PEC). Asserted: it reproduces the record's ``"pec"`` column
at the same band and misses the Cu record by > 1e-5 relative per class.
Predicted: the |cu - pec|/|pec| gap is 6.42e-4 / 2.81e-4 / 4.86e-4
(self / adjacent / opposite). Every band is imported; none is set here.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import FourPortBirdcageGeometry, LeontovichWall, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.lumped import run_lumped_sheet_port_case
from tests.complex_mode import complex_only
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
)
from tests.validation.test_port_birdcage_lumped_column import (
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KEY = "64"
# The negative control's floor (pre-registered, §9 item 324): the PEC Study
# must miss the Cu record by more than this per class.
NEGATIVE_CONTROL_FLOOR = 1e-5


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _rel_means(got, want):
    return {c: abs(got[c] - want[c]) / abs(want[c]) for c in CLASSES}


def _geometry():
    return FourPortBirdcageGeometry(as_hole=True, strip_fraction=study_mod.GATED_WIDTH_FRACTION)


def _identity(res, omega, comm):
    e_res, fields = run_lumped_sheet_port_case(
        res.problem,
        res.port_definitions,
        res.specs,
        facet_tags=res.facet_tags,
        driven_port_id=th14.DRIVEN,
        verbose=False,
        return_fields=True,
        extra_bilinear_terms=res.conductor_terms,
        degree=res.degree,
    )
    comm.Barrier()
    e = fields.e_complex
    sigma_f = fields.sigma_field
    msh = res.mesh
    dx = ufl.Measure("dx", domain=msh, subdomain_data=res.cell_tags)

    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    acct = {"mesh": msh, "facet_tags": res.facet_tags}
    sheets = {sp.port_id: sp.sheet(driven=(sp.port_id == th14.DRIVEN)) for sp in res.specs}
    p_sheets = float(sum(_sheet_field_dissipation_w(acct, sh, e, omega) for sh in sheets.values()))
    p_src = float(_source_power_w(acct, sheets[th14.DRIVEN], e, omega))
    p_total = _vol(dx)
    p_phantom = _vol(dx(PHANTOM_CELL_TAG))
    z_cu = LeontovichWall(study_mod.COPPER_SIGMA).surface_impedance_ohm(omega)
    p_surf = th14._surface_loss_w(msh, res.conductor_facet_tags, z_cu, e, comm)
    resid = abs(p_src - p_sheets - p_phantom - (p_total - p_phantom) - p_surf) / abs(p_src)
    return {
        "p_src": p_src,
        "p_sheets": p_sheets,
        "p_phantom": p_phantom,
        "p_non_phantom": p_total - p_phantom,
        "p_surf": p_surf,
        "residual": resid,
    }


@pytest.fixture(scope="module")
def copper():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"(R) is asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    f_hz = th14.FREQUENCIES_HZ[KEY]
    assert f_hz == FREQUENCY_64_HZ
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    out = {}
    for col, cond in (("cu", LeontovichWall(sigma_s_per_m=study_mod.COPPER_SIGMA)), ("pec", None)):
        study = Study(
            _geometry(), phantom, LumpedPortDrive(), frequency_hz=f_hz, degree=1, conductor=cond
        )
        t0 = time.perf_counter()
        res = study.solve(comm)
        t_solve = time.perf_counter() - t0
        out[col] = {
            "res": res,
            "rec": th14._sweep_record(res.sweep),
            "means": _class_means(res.s_matrix),
            "t": t_solve,
            "omega": study.omega,
        }
        _p0(
            f"\n[WF-16 s7] {col}: cells {res.cells}, mesh {res.mesh_time:.1f} s, sweep "
            f"{res.sweep_time:.1f} s, total {t_solve:.1f} s, -n {comm.size}"
        )
        _p0(
            f"    {col} class means (10 digits): "
            + ", ".join(
                f"{c} {out[col]['means'][c].real:+.9e}{out[col]['means'][c].imag:+.9e}j"
                for c in CLASSES
            )
        )
    t0 = time.perf_counter()
    out["power"] = _identity(out["cu"]["res"], out["cu"]["omega"], comm)
    _p0(
        f"[WF-16 s7] identity field solve {time.perf_counter() - t0:.1f} s; power "
        f"{out['power']}"
    )
    for col in ("cu", "pec"):
        rel = _rel_means(out[col]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY][col])
        _p0(
            f"[WF-16 s7] (R) {col} vs record, relative: "
            + "  ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
            + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
        )
    cu_m, pec_m = out["cu"]["means"], out["pec"]["means"]
    _p0(
        "[WF-16 s7] cu vs pec class-mean relative gap |cu-pec|/|pec| (fraction, PREDICTED "
        "6.42e-4 / 2.81e-4 / 4.86e-4): "
        + "  ".join(f"{c} {abs(cu_m[c] - pec_m[c]) / abs(pec_m[c]):.3e}" for c in CLASSES)
    )
    return out


def test_src_constants_equal_test_values():
    """Ruling 324-4: the src constants are the test values, by import."""
    assert study_mod.GATED_WIDTH_FRACTION == GATED_WIDTH_FRACTION
    assert study_mod.COPPER_SIGMA == th14.COPPER_SIGMA
    assert study_mod.OUTER_BOX_TAG == th14.OUTER_BOX_TAG
    assert LeontovichWall().sigma_s_per_m == th14.COPPER_SIGMA


def test_leontovich_needs_hole_route():
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    study = Study(
        FourPortBirdcageGeometry(),
        phantom,
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        conductor=LeontovichWall(),
    )
    with pytest.raises(ValueError, match="as_hole"):
        study_mod._conductor_setup(study, {}, MPI.COMM_WORLD)


@complex_only
def test_hole_cell_record(copper):
    for col in ("cu", "pec"):
        assert copper[col]["res"].cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD, col


@complex_only
@pytest.mark.parametrize("col", ["cu", "pec"])
def test_r_class_means_reproduce_record(copper, col):
    rel = _rel_means(copper[col]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY][col])
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (col, c, rel[c])


@complex_only
def test_negative_control_pec_misses_cu_record(copper):
    rel = _rel_means(copper["pec"]["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    for c in CLASSES:
        assert rel[c] > NEGATIVE_CONTROL_FLOOR, (c, rel[c])


@complex_only
def test_b_copper_surface_loss_identity(copper):
    p = copper["power"]
    assert p["p_surf"] > 0.0, p
    assert p["residual"] <= DISCRETE_IDENTITY_RTOL, p


@complex_only
@pytest.mark.parametrize("col", ["cu", "pec"])
def test_g_gates(copper, col):
    sw = copper[col]["rec"]
    assert sw["reciprocity"] <= RECIPROCITY_BAND, (col, sw["reciprocity"])
    assert float(np.max(sw["sigma"])) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (col, sw["sigma"])
    for name, value in sw["spreads"].items():
        assert value <= ADJACENT_SPREAD_BAND, (col, name, value)
