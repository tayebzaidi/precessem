"""Example (`EX-108`): the ``Study`` on a copper birdcage -- the F-small four-port coil
with a Leontovich wall on the coil hole, against the PEC control.

`WF-16` step 7 added ``LeontovichWall`` and ``FourPortBirdcageGeometry(as_hole=True,
strip_fraction=...)`` to the ``Study`` API.  This example runs the same F-small
birdcage on exact port rectangles (kappa off) twice -- once on copper, once with
``conductor=None`` (the hole wall left natural PEC) -- at degree 1, 64 MHz, ``-n 2``,
and reads the copper ``Study``'s loss partition (sheets / phantom / copper) and where
on the coil the copper loss sits.

Every record, band and helper is imported from the gate modules (the ANS-1 rule); the
example states no sigma, epsilon_r, frequency, band, power or cell count of its own.

Asserted:

* (W) ``comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`` before any comparison;
* (M) both Studies mesh to ``ANS6_EXACT_STRIP_HOLE_CELL_RECORD`` cells;
* (R) the copper and PEC class means equal `ANS-6` step 4's registered means at
  ``LEG_D0_REPRODUCTION_BAND``;
* (b) the copper power identity: ``p_surf > 0`` and residual <= ``DISCRETE_IDENTITY_RTOL``;
* (g) `TH-14`'s gates (reciprocity, passivity, class spread) on both columns;
* (L) the per-cell copper loss field sums to ``p_surf`` within ``DISCRETE_IDENTITY_RTOL``.

Negative control (asserted): (c1) the PEC ``Study`` cannot stand in for copper -- its
class means miss the copper record by more than ``NEGATIVE_CONTROL_FLOOR`` per class.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:16 -n 2 -t 300

Not claimed: an absolute S11 / Z_in, a resonance, tuning or Q figure; degree 2;
10 / 128 MHz; a surface loss density; an accuracy statement on the copper loss.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem, la

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import (  # noqa: E402
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial  # noqa: E402
from fem_em_solver.ports.lumped import run_lumped_sheet_port_case  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.study import COPPER_SIGMA, GATED_WIDTH_FRACTION  # noqa: E402

from tests.package.test_wf16_study_copper import (  # noqa: E402
    KEY,
    NEGATIVE_CONTROL_FLOOR,
    _identity,
    _rel_means,
)
from tests.validation import test_th14_birdcage_copper as th14  # noqa: E402
from tests.validation.test_ans6_exact_strip import (  # noqa: E402
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_power_identity import (  # noqa: E402
    DISCRETE_IDENTITY_RTOL,
)
from tests.validation.test_lossy_sphere_fullwave import (  # noqa: E402
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
)
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_lumped_narrowed_sheet import (  # noqa: E402
    GATED_WIDTH_FRACTION as GATE_WIDTH_FRACTION,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_16_study_copper_birdcage"
DEGREE = 1  # the record's degree (ANS-6 step 4 P1); the Study default is 2.


def _coil_surface_loss_per_cell_w(msh, facet_tags, z_s, e, comm):
    """Facet-integrated copper loss ``0.5 Re(1/Z_s)|n x E|^2`` per *adjacent cell*, W.

    Not a density.  The same integrand ``th14._surface_loss_w`` reduces to a scalar,
    with a DG0 test function inserted so each cell touching facet tag 401 keeps its
    own share (the ``ports:18`` construction, rewritten here: examples never import
    each other).  A DG0 test function on an exterior facet takes its one adjacent
    cell's value, so no interior-restriction ambiguity arises.
    """
    v0 = fem.functionspace(msh, ("DG", 0))
    v = ufl.TestFunction(v0)
    ds_c = ufl.Measure(
        "ds", domain=msh, subdomain_data=facet_tags, subdomain_id=BIRDCAGE_CONDUCTOR_SURFACE_TAG
    )
    n = ufl.FacetNormal(msh)
    re_y = fem.Constant(msh, default_scalar_type(0.5 * np.real(1.0 / complex(z_s))))
    # conj(v): complex-mode arity needs the test function conjugated; v is real DG0.
    form = fem.form(re_y * ufl.inner(ufl.cross(n, e), ufl.cross(n, e)) * ufl.conj(v) * ds_c)
    b = fem.assemble_vector(form)
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    loss = fem.Function(v0, name="coil_surface_loss_per_cell_W")
    loss.x.array[:] = np.real(b.array)
    loss.x.scatter_forward()
    size_local = v0.dofmap.index_map.size_local
    total = float(comm.allreduce(np.sum(np.real(b.array[:size_local])), op=MPI.SUM))
    return loss, total


def _means_line(means):
    return ", ".join(f"{c} {means[c].real:+.9e}{means[c].imag:+.9e}j" for c in CLASSES)


def _rel_line(rel):
    return "  ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- (W) the record holds at its own width -----------------------------
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )
    # The package constants are the gate modules' values (ruling 324-4), by import.
    assert GATED_WIDTH_FRACTION == GATE_WIDTH_FRACTION
    assert COPPER_SIGMA == th14.COPPER_SIGMA
    assert th14.FREQUENCIES_HZ[KEY] == FREQUENCY_64_HZ

    # ---- the two Studies (a fresh geometry each) ---------------------------
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    geometry = FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)
    study_cu = Study(
        geometry, phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ,
        degree=DEGREE, conductor=LeontovichWall(COPPER_SIGMA),
    )
    t0 = time.perf_counter()
    cu = study_cu.solve(comm)
    t_cu = time.perf_counter() - t0

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    msh = cu.mesh
    local_tags = np.unique(cu.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    air_tags = tuple(t for t in all_tags if t != PHANTOM_CELL_TAG)
    region_names = {PHANTOM_CELL_TAG: "gel phantom"}
    for t in air_tags:
        region_names[t] = "air (outer box PEC)"
    write_setup_figure(
        msh,
        cu.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=air_tags,
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title=(
            "mri:16 -- Study on the copper F-small birdcage: Leontovich wall on the "
            "coil hole (WF-16 step 7, degree 1)"
        ),
        comm=comm,
    )

    geometry_pec = FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)
    study_pec = Study(
        geometry_pec, phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ,
        degree=DEGREE, conductor=None,
    )
    t0 = time.perf_counter()
    pec = study_pec.solve(comm)
    t_pec = time.perf_counter() - t0

    # ---- the copper power partition (the gate's own computation) -----------
    t0 = time.perf_counter()
    power = _identity(cu, study_cu.omega, comm)
    t_id = time.perf_counter() - t0

    # ---- one more P1 field solve for the XDMF (`_identity` returns no field)
    t0 = time.perf_counter()
    _e_res, fields = run_lumped_sheet_port_case(
        cu.problem, cu.port_definitions, cu.specs, facet_tags=cu.facet_tags,
        driven_port_id=th14.DRIVEN, verbose=False, return_fields=True,
        extra_bilinear_terms=cu.conductor_terms, degree=cu.degree,
    )
    t_field = time.perf_counter() - t0
    e = fields.e_complex

    z_s = LeontovichWall(COPPER_SIGMA).surface_impedance_ohm(study_cu.omega)
    loss_field, loss_total = _coil_surface_loss_per_cell_w(
        msh, cu.conductor_facet_tags, z_s, e, comm
    )
    # Every collective on every rank: |E| into P1 (the mesh is affine, so the N1curl
    # Function is interpolated, not written).
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()

    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cu.cell_tags,
        {"E_abs": e_abs, "coil_surface_loss_per_cell_W": loss_field},
        comm=comm,
        facet_tags=cu.conductor_facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- readings ----------------------------------------------------------
    means_cu = _class_means(cu.s_matrix)
    means_pec = _class_means(pec.s_matrix)
    rec = ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]
    rel_cu = _rel_means(means_cu, rec["cu"])
    rel_pec = _rel_means(means_pec, rec["pec"])
    c1 = _rel_means(means_pec, rec["cu"])  # the asserted form (relative to the Cu record)
    gap_gate = {c: abs(means_cu[c] - means_pec[c]) / abs(means_pec[c]) for c in CLASSES}
    gap_cu = {c: abs(means_cu[c] - means_pec[c]) / abs(means_cu[c]) for c in CLASSES}
    sweeps = {"cu": th14._sweep_record(cu.sweep), "pec": th14._sweep_record(pec.sweep)}

    p_in = power["p_src"] - power["p_sheets"]
    frac_coil = power["p_surf"] / p_in
    frac_phantom = power["p_phantom"] / p_in
    loss_rel = abs(loss_total - power["p_surf"]) / abs(power["p_surf"])

    # ---- prints (rank 0, no collective inside) -----------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-108 -- the Study on a copper birdcage (WF-16 step 7)")
        print("=" * 72)
        print(
            f"[study] degree {study_cu.degree}, -n {comm.size}, f = "
            f"{study_cu.frequency_hz:.6e} Hz, copper sigma = {COPPER_SIGMA:g} S/m, "
            f"Z_s = {z_s:.6e} ohm"
        )
        print(
            f"[mesh] copper {cu.cells} cells, PEC {pec.cells} cells (record "
            f"{ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); copper Study {t_cu:.2f} s (mesh "
            f"{cu.mesh_time:.2f} s, sweep {cu.sweep_time:.2f} s), PEC Study {t_pec:.2f} s, "
            f"identity solve {t_id:.2f} s, field solve {t_field:.2f} s"
        )
        print(f"[(W)] comm.size {comm.size} == record width {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}")
        print(f"[cu] class means (10 digits): {_means_line(means_cu)}")
        print(f"[pec] class means (10 digits): {_means_line(means_pec)}")
        print(
            f"[(R)] cu vs record, relative: {_rel_line(rel_cu)}; pec vs record: "
            f"{_rel_line(rel_pec)} (band {LEG_D0_REPRODUCTION_BAND:.0e})"
        )
        print(
            "[gap] |cu-pec|/|pec| (the gate's form): "
            + "  ".join(f"{c} {gap_gate[c]:.3e}" for c in CLASSES)
            + "; |cu-pec|/|cu| (relative to the Cu mean): "
            + "  ".join(f"{c} {gap_cu[c]:.3e}" for c in CLASSES)
        )
        print(
            f"[(c1)] PEC Study vs the Cu record, relative (asserted form): "
            f"{_rel_line(c1)} (floor {NEGATIVE_CONTROL_FLOOR:g}, "
            f"{min(c1.values()) / NEGATIVE_CONTROL_FLOOR:.1f}x; "
            f"{min(c1.values()) / LEG_D0_REPRODUCTION_BAND:.1e}x the record band)"
        )
        for col in ("cu", "pec"):
            sw = sweeps[col]
            print(
                f"[(g)] {col}: reciprocity {sw['reciprocity']:.3e} (band {RECIPROCITY_BAND:g}); "
                f"sigma_max {float(np.max(sw['sigma'])):.12f} (<= 1 + "
                f"{PASSIVITY_SIGMA_TOLERANCE:g}); spreads "
                + ", ".join(f"{k} {v * 100:.4f}%" for k, v in sw["spreads"].items())
                + f" (band {ADJACENT_SPREAD_BAND:g})"
            )
        print(
            f"[(b)] copper power identity: p_surf {power['p_surf']!r} W, residual "
            f"{power['residual']:.3e} (DISCRETE_IDENTITY_RTOL {DISCRETE_IDENTITY_RTOL:g})"
        )
        print(
            "[partition] p_src {p_src:.9e} W, p_sheets {p_sheets:.9e} W, p_phantom "
            "{p_phantom:.9e} W, p_non_phantom {p_non_phantom:.3e} W, p_surf {p_surf:.9e} W"
            .format(**power)
        )
        print(
            f"[partition] P_in = p_src - p_sheets = {p_in:.6e} W; P_coil/P_in "
            f"{frac_coil:.7f}; P_phantom/P_in {frac_phantom:.7f}"
        )
        print(
            f"[(L)] sum of coil_surface_loss_per_cell_W over owned dofs {loss_total!r} W "
            f"(facet-integrated loss per adjacent cell, not a density); vs p_surf relative "
            f"{loss_rel:.3e} (DISCRETE_IDENTITY_RTOL {DISCRETE_IDENTITY_RTOL:g})"
        )
        print(f"[paraview] wrote {xdmf_path}")

    # ---- asserts -----------------------------------------------------------
    assert cu.cells == pec.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD, (
        f"(M) cells cu {cu.cells} pec {pec.cells} record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}"
    )
    for c in CLASSES:
        assert rel_cu[c] <= LEG_D0_REPRODUCTION_BAND, f"(R) cu {c} {rel_cu[c]:.3e}"
        assert rel_pec[c] <= LEG_D0_REPRODUCTION_BAND, f"(R) pec {c} {rel_pec[c]:.3e}"
        assert c1[c] > NEGATIVE_CONTROL_FLOOR, f"(c1) blind: {c} {c1[c]:.3e}"
    assert power["p_surf"] > 0.0, f"(b) {power}"
    assert power["residual"] <= DISCRETE_IDENTITY_RTOL, f"(b) {power}"
    for col, sw in sweeps.items():
        assert sw["reciprocity"] <= RECIPROCITY_BAND, f"(g) {col} reciprocity {sw['reciprocity']:.3e}"
        assert float(np.max(sw["sigma"])) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (
            f"(g) {col} sigma {sw['sigma']}"
        )
        for name, value in sw["spreads"].items():
            assert value <= ADJACENT_SPREAD_BAND, f"(g) {col} spread {name} {value:.3e}"
    assert loss_rel <= DISCRETE_IDENTITY_RTOL, f"(L) per-cell sum vs p_surf {loss_rel:.3e}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
