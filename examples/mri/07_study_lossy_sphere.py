"""Example (`EX-91`): the first study through the installed API -- a ``Study`` on the
lossy sphere reproduces `MAT-4`'s SAR.

`WF-15`'s ``Study`` composes a geometry generator, a gelled-saline material, an
imposed uniform field, a solve and a SAR result. This example builds it on
`MAT-4` step 1's lossy sphere at the **default degree** (no ``degree=`` keyword,
so 2) and reads its mean SAR against `WF-15`'s degree-2 record and against
`MAT-4`'s closed form.

Every constant is imported from the gate modules (the ANS-1 rule); the example
states no sigma, epsilon_r, frequency, field or radius of its own -- it reads
them from the built ``Study``.

Asserted:

* (Deg) ``study.degree == 2`` and the solved N1curl space is degree 2;
* (R)   |SAR - WF15_DEG2_RECORD| / record <= ``L0_BAND``;
* (G)   |SAR - closed form| / closed form <= ``MAT4_BAND``.

Negative controls (both asserted):

* (c1) the sigma-shifted ``Study`` (``material="mid"``, degree 1) reproduces
  `MAT-1`'s ``.4e`` print and misses the sigma = 0.57 closed form by more than
  ``MAT4_BAND``;
* (c2) the (R) comparator separates degrees: ``MAT4_DEG1_RECORD`` against
  ``WF15_DEG2_RECORD`` differs by more than ``L0_BAND`` (two imported constants,
  no solve).

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:7 -n 2 -t 590

Not claimed: anything `WF-15` does not claim -- no coil, no convergence or
h-statement, no local or 10 g SAR accuracy, no degree-1 production use.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import Study  # noqa: E402,F401  (the API under test, via _study)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.materials.test_phantom_material_model import MAT4_BAND  # noqa: E402
from tests.package.test_wf15_study_lossy_sphere import (  # noqa: E402
    GEOMETRY,
    MAT1_MID_PRINT,
    MAT4_DEG1_RECORD,
    _closed,
    _study,
)
from tests.package.test_wf17_first_study_tutorial import WF15_DEG2_RECORD  # noqa: E402
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_07_study_lossy_sphere"
AIR_TAG = 2


def _global_dofs(V) -> int:
    return int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- the study at the default degree -----------------------------------
    study = _study()  # no degree= keyword
    t0 = time.perf_counter()
    res = study.solve(comm)
    solve_seconds = time.perf_counter() - t0
    sar_dict = study.mean_sar()
    sar = float(np.real(sar_dict["mean_sar_w_per_kg"]))

    msh = res.mesh
    sigma_value = study.material.sigma
    eps_value = study.material.epsilon_r
    freq = study.frequency_hz

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        res.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            study.geometry.sphere_tag: f"gel sphere (sigma {sigma_value:g} S/m)",
            AIR_TAG: "air box (imposed-field walls)",
        },
        hide_tags=(AIR_TAG,),
        translucent_tags=(study.geometry.sphere_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:7 -- Study on the lossy sphere (WF-15, default degree)",
        comm=comm,
    )

    V = res.fields.e_complex.function_space
    element_degree = V.ufl_element().degree
    element_degree = element_degree() if callable(element_degree) else element_degree
    n_cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    ndofs = _global_dofs(V)

    closed = _closed(sigma_value, eps_value)["sar_w_per_kg"]
    rel_record = abs(sar - WF15_DEG2_RECORD) / WF15_DEG2_RECORD
    rel_closed = abs(sar - closed) / closed

    # ---- (c1): the sigma-shifted Study, degree 1 ---------------------------
    study_mid = _study(material="mid", degree=1)
    study_mid.solve(comm)
    sar_mid = float(np.real(study_mid.mean_sar()["mean_sar_w_per_kg"]))
    mid_print = f"{sar_mid:.4e}"
    rel_mid_closed = abs(sar_mid - closed) / closed

    # ---- (c2): the (R) comparator separates degrees (no solve) -------------
    c2 = abs(MAT4_DEG1_RECORD / WF15_DEG2_RECORD - 1.0)

    # ---- output fields (every collective on every rank) --------------------
    e = res.fields.e_complex
    p1 = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()

    sigma_field = res.fields.sigma_field
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar_field = fem.Function(dg0, name="SAR")
    sar_field.interpolate(
        fem.Expression(
            sigma_field * ufl.real(ufl.inner(e, e)) / (2.0 * study.density_kg_m3),
            dg0.element.interpolation_points,
        )
    )
    sar_field.x.array[:] = np.real(sar_field.x.array)
    sar_field.x.scatter_forward()
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()

    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        res.cell_tags,
        {"E_abs": e_abs, "SAR": sar_field, "sigma": sigma_out},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) -----------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-91 -- a Study on the lossy sphere at the default degree (WF-15)")
        print("=" * 72)
        print(
            f"[study] sigma = {sigma_value} S/m, eps_r = {eps_value}, f = {freq:.6e} Hz; "
            f"Study.degree = {study.degree}, result.degree = {res.degree}, "
            f"space degree {element_degree}"
        )
        print(
            f"[mesh] {n_cells} cells, {ndofs} N1curl dofs (printed; version-tagged "
            f"74020 / 482158); Study solve+mesh {solve_seconds:.1f} s"
        )
        print(f"[(Deg)] study.degree={study.degree} element degree={element_degree}")
        print(
            f"[(R)] SAR = {sar!r} ({sar:.10e}); record {WF15_DEG2_RECORD:.10e}; "
            f"rel {rel_record:.3e} (band L0_BAND {L0_BAND:g})"
        )
        print(
            f"[(G)] closed form {closed:.6e}; rel {rel_closed:.4%} (band MAT4_BAND {MAT4_BAND:g})"
        )
        print(
            f"[(c1)] mid preset (sigma {study_mid.material.sigma}, degree {study_mid.degree}): "
            f"SAR = {sar_mid:.10e} ('{mid_print}' vs MAT-1 '{MAT1_MID_PRINT}'); "
            f"rel to closed(0.57) {rel_mid_closed:.4%} (must exceed MAT4_BAND {MAT4_BAND:g})"
        )
        print(
            f"[(c2)] |MAT4_DEG1_RECORD/WF15_DEG2_RECORD - 1| = {c2:.4e} "
            f"(must exceed L0_BAND {L0_BAND:g}; degree-1 SAR vs this record measured "
            f"3.494e-02 in 20261006T004414Z_WF-17-W1-planted-deg1.log:161)"
        )
        print(f"[paraview] wrote {xdmf_path}")

    assert study.degree == 2 and res.degree == 2, "(Deg) default degree is not 2"
    assert element_degree == 2, f"(Deg) solved space degree {element_degree}"
    assert rel_record <= L0_BAND, f"(R) rel {rel_record:.3e} > {L0_BAND:g}"
    assert rel_closed <= MAT4_BAND, f"(G) rel {rel_closed:.4%} > {MAT4_BAND:g}"
    assert mid_print == MAT1_MID_PRINT, f"(c1) {mid_print} != {MAT1_MID_PRINT}"
    assert rel_mid_closed > MAT4_BAND, f"(c1) blind: {rel_mid_closed:.4%} <= {MAT4_BAND:g}"
    assert c2 > L0_BAND, f"(c2) the (R) comparator does not separate degrees: {c2:.3e}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
