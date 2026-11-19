"""`WF-16` step 11 (§9 item 388, lane 1 defect fix): the field-keeping port
solve (``Study._solve_port_fields``, behind ``Study.sar`` / ``b1_plus`` /
``b1_plus_at``) applies the ``Study``'s conductor exactly as ``Study.solve()``.

Pre-fix the field route (F1) raised ``material_map ... [1]`` on the
``as_hole=True`` copper Study and (F2) silently returned the sigma = 800 S/m
solid fields on ``as_hole=False`` with ``LeontovichWall()``.

- (K) no mesh: the solid-route ``LeontovichWall`` and a non-``LeontovichWall``
  conductor are refused by ``sar`` / ``b1_plus`` as ``solve()`` refuses them;
- (K1) the hole cell record; the Cu result carries the conductor term / facet
  tags, the PEC (``conductor=None``) result carries none;
- (S=) field route S = ``solve()`` S on the same Cu Study at
  ``LEG_D0_REPRODUCTION_BAND``;
- (R) the field route reproduces ``ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]``
  ``cu`` / ``pec`` at ``LEG_D0_REPRODUCTION_BAND``;
- (P) the kept Cu P1 field closes the exact power identity with P_surf > 0 at
  ``DISCRETE_IDENTITY_RTOL``.

Negative controls (rule (e), asserted, backed): (N1) each column misses the
other's record by > ``NEGATIVE_CONTROL_FLOOR`` per class; (N2) the P_surf-blind
residual is >= ``BLIND_FACTOR`` x ``DISCRETE_IDENTITY_RTOL``.  SAR and
``|B1+|`` on copper are printed records only, never asserted against a value.
Every band is imported; none is set here.
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
from tests.complex_mode import complex_only
from tests.package.test_wf16_study_copper import NEGATIVE_CONTROL_FLOOR, _geometry
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
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR

KEY = "64"
CENTRE = (0.0, 0.0, 0.0)


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _rel_means(got, want):
    return {c: abs(got[c] - want[c]) / abs(want[c]) for c in CLASSES}


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)


# Copied, with provenance, from tests/package/test_wf16_study_copper.py:84-110
# (`_identity`, main at e122cead), minus its `run_lumped_sheet_port_case` call:
# the field is the sweep's kept one (the field `sar` / `b1_plus` read).
def _kept_identity(res, omega, comm):
    fields = res.fields[th14.DRIVEN]
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
    blind = abs(p_src - p_sheets - p_total) / abs(p_src)
    return {
        "p_src": p_src,
        "p_sheets": p_sheets,
        "p_phantom": p_phantom,
        "p_non_phantom": p_total - p_phantom,
        "p_surf": p_surf,
        "residual": resid,
        "blind_residual": blind,
    }


@pytest.fixture(scope="module")
def fields():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"(R) is asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    assert th14.FREQUENCIES_HZ[KEY] == FREQUENCY_64_HZ
    out = {}
    # (1) Cu Study A, field route first (pre-fix this raised `material_map ... [1]`).
    a = Study(
        _geometry(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        conductor=LeontovichWall(study_mod.COPPER_SIGMA),
    )
    t0 = time.perf_counter()
    res_f = a._solve_port_fields(comm)
    out["t_f"] = time.perf_counter() - t0
    assert a.result is None
    # (2) solve() on the same Study.
    t0 = time.perf_counter()
    res_s = a.solve(comm)
    out["t_s"] = time.perf_counter() - t0
    # (3) PEC-hole Study B.
    b = Study(
        _geometry(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        conductor=None,
    )
    t0 = time.perf_counter()
    res_p = b._solve_port_fields(comm)
    out["t_p"] = time.perf_counter() - t0
    _p0(
        f"\n[WF-16 s11] cells cu {res_f.cells} / pec {res_p.cells}; field route cu "
        f"{out['t_f']:.1f} s (mesh {res_f.mesh_time:.1f}, sweep {res_f.sweep_time:.1f}), "
        f"solve() cu {out['t_s']:.1f} s, field route pec {out['t_p']:.1f} s "
        f"(mesh {res_p.mesh_time:.1f}, sweep {res_p.sweep_time:.1f}), -n {comm.size}"
    )
    # (4) the kept-field identity.
    t0 = time.perf_counter()
    power = _kept_identity(res_f, a.omega, comm)
    _p0(
        f"[WF-16 s11] kept-field identity {time.perf_counter() - t0:.1f} s: "
        f"P_src {power['p_src']:.9e} W, P_sheets {power['p_sheets']:.9e} W, "
        f"P_phantom {power['p_phantom']:.9e} W, P_non_phantom {power['p_non_phantom']:.9e} W, "
        f"P_surf {power['p_surf']:.9e} W, residual {power['residual']:.3e} "
        f"(band {DISCRETE_IDENTITY_RTOL:.0e}), blind residual {power['blind_residual']:.4e}"
    )
    _p0(
        "[WF-16 s11] direct-route prior (item360-W9-copper.log:1823): residual 8.507e-13, "
        "P_src 4.806472101e-03 W, P_surf 1.516873149e-06 W, P_phantom 1.870050049e-06 W"
    )
    # (5) records (collective on every rank).
    t0 = time.perf_counter()
    rec = {}
    for col, st in (("cu", a), ("pec", b)):
        sar = st.sar(mass_g=10.0, centres=[CENTRE], drives=["P1"], comm=comm)
        b1 = st.b1_plus_at([CENTRE], drives=["P1"], comm=comm)
        rec[col] = {"sar": sar, "b1": b1}
        _p0(
            f"[WF-16 s11] record {col}: SAR10g(P1, centre) "
            f"{sar['P1'][0]['averaged_sar_w_per_kg']:.9e} W/kg (mass "
            f"{sar['P1'][0]['mass_kg']:.6e} kg); |B1+|(P1, centre) "
            f"{float(b1['P1'][0][0]):.9e} T valid {bool(b1['P1'][1][0])} (per 1 V, printed only)"
        )
    _p0(f"[WF-16 s11] records {time.perf_counter() - t0:.1f} s")
    out.update(
        a=a,
        b=b,
        res_f=res_f,
        res_s=res_s,
        res_p=res_p,
        power=power,
        rec=rec,
        means_f=_class_means(res_f.s_matrix),
        means_p=_class_means(res_p.s_matrix),
    )
    return out


def test_wf16_s11_conductor_refused_before_mesh(monkeypatch):
    """(K) the field route refuses where `solve()` refuses, before any mesh."""
    monkeypatch.setattr(FourPortBirdcageGeometry, "build", lambda self, comm: {"mesh": None})
    solid = Study(
        FourPortBirdcageGeometry(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        conductor=LeontovichWall(),
    )
    with pytest.raises(ValueError, match="as_hole"):
        solid.sar(centres=[CENTRE])
    assert getattr(solid, "fields_result", None) is None
    with pytest.raises(ValueError, match="as_hole"):
        solid.b1_plus()
    assert getattr(solid, "fields_result", None) is None
    with pytest.raises(ValueError, match="as_hole"):
        solid.solve()
    other = Study(
        FourPortBirdcageGeometry(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        conductor=object(),
    )
    with pytest.raises(TypeError):
        other.sar(centres=[CENTRE])
    assert getattr(other, "fields_result", None) is None


@complex_only
def test_wf16_s11_cells_and_terms(fields):
    res_f, res_p = fields["res_f"], fields["res_p"]
    assert res_f.cells == res_p.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert res_f.conductor_terms is not None and len(res_f.conductor_terms) == 1
    assert res_f.conductor_facet_tags is not None
    assert res_p.conductor_terms is None
    assert res_p.conductor_facet_tags is None


@complex_only
def test_wf16_s11_field_route_equals_solve(fields):
    s_f = np.asarray(fields["res_f"].s_matrix)
    s_s = np.asarray(fields["res_s"].s_matrix)
    rel = float(np.linalg.norm(s_f - s_s) / np.linalg.norm(s_s))
    _p0(
        f"\n[WF-16 s11] (S=) ||S_field - S_solve||_F / ||S_solve||_F {rel:.3e} "
        f"(band {LEG_D0_REPRODUCTION_BAND:.0e}; predicted <= 1e-14)"
    )
    assert rel <= LEG_D0_REPRODUCTION_BAND, rel


@complex_only
@pytest.mark.parametrize("col", ["cu", "pec"])
def test_wf16_s11_records(fields, col):
    """(R) the field route reproduces the fixture's record; SAR / |B1+| printed."""
    means = fields["means_f"] if col == "cu" else fields["means_p"]
    rel = _rel_means(means, ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY][col])
    _p0(
        f"\n[WF-16 s11] (R) {col} field route vs record, relative: "
        + "  ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (col, c, rel[c])
    rec = fields["rec"][col]
    assert list(rec["sar"]) == ["P1"] and len(rec["sar"]["P1"]) == 1
    assert bool(rec["b1"]["P1"][1][0])


@complex_only
def test_wf16_s11_kept_field_power_identity(fields):
    p = fields["power"]
    assert p["p_surf"] > 0.0, p
    assert p["residual"] <= DISCRETE_IDENTITY_RTOL, p


@complex_only
def test_wf16_s11_negative_wrong_column(fields):
    """(N1) asserted, backed: each column misses the other's record."""
    for col, means, other in (("cu", fields["means_f"], "pec"), ("pec", fields["means_p"], "cu")):
        rel = _rel_means(means, ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY][other])
        _p0(
            f"\n[WF-16 s11] (N1) {col} vs {other} record: "
            + "  ".join(
                f"{c} {rel[c]:.3e} ({rel[c] / NEGATIVE_CONTROL_FLOOR:.1f}x floor)" for c in CLASSES
            )
        )
        for c in CLASSES:
            assert rel[c] > NEGATIVE_CONTROL_FLOOR, (col, other, c, rel[c])


@complex_only
def test_wf16_s11_negative_surface_blind(fields):
    """(N2) asserted, backed: the P_surf-blind residual clears the 10x bar."""
    p = fields["power"]
    bar = BLIND_FACTOR * DISCRETE_IDENTITY_RTOL
    _p0(
        f"\n[WF-16 s11] (N2) blind residual {p['blind_residual']:.4e} = "
        f"{p['blind_residual'] / DISCRETE_IDENTITY_RTOL:.1f}x DISCRETE_IDENTITY_RTOL "
        f"(bar {bar:.0e}; ceiling ≈ 3.1559e-04)"
    )
    assert p["blind_residual"] >= bar, p
