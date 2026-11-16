"""Example (`EX-72`): a metal sphere in the gel -- the dipole it adds and the
SAR it concentrates.

`IMP-1` step 1 gates a perfectly conducting sphere (a = 5 mm) cut out of a
10 mm half-width box as a hole, in `MAT-4`'s gelled-saline host (sigma =
0.57 S/m, eps_r = 78) at 64 MHz, degree 2, under an imposed quasi-static
field E0 zhat. `th:9` shows the same hole lossless at 4.77 MHz and reads only
beta; this example adds the lossy host and the quantity an implant is about:
the shell-mean SAR enhancement over a < r < 1.5 a,

    Phi = <|E|^2>/E0^2 = 1 + 2 beta^2 (a/c)^3 = 1.5926   (c = 1.5 a, beta = 1),

and the pole / equator SAR pattern at r = 1.2 a.

**It asserts, it does not merely render.** Every record, band and helper is
*imported* from the module that gates the capability
(``tests/validation/test_imp1_pec_sphere_lossy_host.py``), never restated:

* |beta - 1| <= `TH15_BETA_BAND`; |Phi_h / Phi_cf,h - 1| <= `MAT4_BAND`;
  cavity dofs present and max |E| on them < 1e-12 -- each printed beside the
  gate's own record.
* Negative control, asserted on the same mesh: the cavity left natural
  (``pec_facet_tags=(OUTER_BOUNDARY_TAG,)``) is a void, not a conductor, and
  misses both bands.

Scope: a/6 rung, degree 2, imposed quasi-static field. No coil, no tissue, no
heating, no ASTM F2182 claim.

Run it through the example runner (the ``imp:`` group sources the complex
build automatically; a real build raises)::

    ./scripts/run_examples.sh -e imp:1 -n 2 -t 180

Output lands in ``examples/implants/paraview_output``: open
``implants_01_pec_sphere_in_gel_combined.xdmf`` and colour by ``sar_ratio``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

from fem_em_solver.io.paraview_utils import (
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure

# The gated fixture lives in the test that gates `IMP-1` step 1; the runner
# puts only ``src`` on PYTHONPATH, so the repo root goes on sys.path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation.test_imp1_pec_sphere_lossy_host import (  # noqa: E402
    BETA_PEC,
    BETA_VOID,
    CAVITY_TAG,
    DEGREE,
    E0,
    FREQUENCY_HZ,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    PHI_CLOSED_FORM_PEC,
    SIGMA_HIGH,
    SPHERE_RADIUS,
    TH15_BETA_BAND,
    _hole_mesh,
    _sar_ratio_records,
    _solve,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_01_pec_sphere_in_gel"


def _paraview_fields(msh, fields):
    """CG1 ``E_real`` / ``E_magnitude`` and DG0 ``sar_ratio`` = |E|^2 / E0^2."""
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_cg = fem.Function(v_cg, name="E_phasor")
    e_cg.interpolate(fields.e_complex)
    e_cg.x.scatter_forward()

    e_re = fem.Function(v_cg, name="E_real")
    e_re.x.array[:] = np.real(e_cg.x.array)
    e_re.x.scatter_forward()

    s_cg = fem.functionspace(msh, ("Lagrange", 1))
    e_mag = fem.Function(s_cg, name="E_magnitude")
    components = np.abs(e_cg.x.array.reshape(-1, 3))
    e_mag.x.array[:] = np.sqrt(np.sum(components * components, axis=1))
    e_mag.x.scatter_forward()

    # local SAR = sigma |E|^2 / (2 rho) in a homogeneous host, so |E|^2/E0^2
    # *is* SAR/SAR0 (SAR0 = sigma E0^2 / (2 rho)).
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar_ratio = fem.Function(dg0, name="sar_ratio")
    e = fields.e_complex
    expr = ufl.real(ufl.inner(e, e)) / (E0**2)
    sar_ratio.interpolate(fem.Expression(expr, dg0.element.interpolation_points))
    sar_ratio.x.scatter_forward()
    return e_re, e_mag, sar_ratio


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    if comm.rank == 0:
        print("=" * 72)
        print("EX-72 -- a metal sphere in the gel: the dipole it adds, the SAR it concentrates")
        print("=" * 72)
        print(
            f"\n[geometry] PEC sphere a = {SPHERE_RADIUS} m cut out as a hole "
            f"(cavity = facet tag {CAVITY_TAG}), host sigma = {SIGMA_HIGH} S/m, "
            f"f = {FREQUENCY_HZ / 1e6:.3f} MHz, degree {DEGREE}, imposed uniform E0"
            f"\n[closed form] beta = {BETA_PEC:.1f}; shell-mean SAR enhancement "
            f"Phi(1.5a) = {PHI_CLOSED_FORM_PEC:.4f}",
            flush=True,
        )

    # One mesh cache: the gate's own `_hole_mesh()` serves the solve, the
    # figure and the XDMF.
    msh, cell_tags, facet_tags = _hole_mesh()

    solve_started = time.perf_counter()
    anchor = _solve((OUTER_BOUNDARY_TAG, CAVITY_TAG), DEGREE)
    sar_records = _sar_ratio_records(anchor["fields"])  # collective: every rank
    solve_seconds = time.perf_counter() - solve_started

    beta = anchor["beta"]
    beta_miss = abs(beta - BETA_PEC)
    phi_h, phi_cf_h = anchor["phi_h"], anchor["phi_cf_h"]
    phi_miss = abs(phi_h / phi_cf_h - 1.0)

    # EX-57 setup figure: cell tag 2 only (the sphere is a hole, there is no
    # tag 1); clip y < 0 opens the block to show the cavity, slice at z = 0.
    # No-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={2: "gel host sigma 0.57 S/m (sphere cut out as a hole)"},
        clip_normal=(0.0, -1.0, 0.0),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:1 -- PEC sphere (a = 5 mm) cut out of a gel host: cavity wall is a facet tag",
        comm=comm,
    )

    if comm.rank == 0:
        print(
            f"\n[solve] {anchor['ncells']} cells (IMP-1 a/6 rung), anchor + SAR "
            f"records in {solve_seconds:.1f} s"
            f"\n[anchor] fitted beta = {beta:.6f}, |beta - 1| = {beta_miss:.4%} "
            f"vs TH15_BETA_BAND = {TH15_BETA_BAND:.4%} "
            f"(gate record: beta 0.985596, `20261001T034025Z_IMP-1-land-W3.log:371`)"
            f"\n[anchor] Phi_h = {phi_h:.6f}, Phi_cf,h = {phi_cf_h:.6f}, "
            f"|Phi_h/Phi_cf,h - 1| = {phi_miss:.4%} vs MAT4_BAND = {MAT4_BAND:.0%} "
            f"(gate record: 1.606456 vs 1.605264); continuous closed form "
            f"{PHI_CLOSED_FORM_PEC:.4f}"
            f"\n[structure] cavity dofs (reduced) = {anchor['n_cavity']}, max |E| "
            f"on them = {anchor['cavity_max']:.3e} (gate record: 5670, 0)",
            flush=True,
        )
        for theta_deg, ratio, closed in sar_records:
            label = "pole" if theta_deg < 45 else "equator"
            record = "4.5288 (cf 4.6421)" if label == "pole" else "0.2072 (0.1898)"
            print(
                f"[SAR] SAR/SAR0 at r = 1.2a, polar {theta_deg:.0f} deg ({label}): "
                f"{ratio:.4f} (quasi-static closed form {closed:.4f}); gate record "
                f"{record}",
                flush=True,
            )

    assert anchor["n_cavity"] > 0, "no dofs located on cavity facet tag 2"
    assert anchor["cavity_max"] < 1e-12, (
        f"max |E| on the cavity dofs is {anchor['cavity_max']:.3e}, not zero"
    )
    assert beta_miss <= TH15_BETA_BAND, (
        f"beta = {beta:.6f} misses 1 by {beta_miss:.4%} > {TH15_BETA_BAND:.4%} "
        "(imported `TH-15` band) -- an example/test divergence, not a tolerance to move"
    )
    assert phi_miss <= MAT4_BAND, (
        f"Phi_h = {phi_h:.6f} vs Phi_cf,h = {phi_cf_h:.6f}: miss {phi_miss:.4%} "
        f"> {MAT4_BAND:.0%} (imported `MAT-4` band) -- an example/test divergence"
    )

    # ---- negative control: asserted, same mesh -------------------------------
    ctl = _solve((OUTER_BOUNDARY_TAG,), DEGREE)
    delta_beta = abs(ctl["beta"] - BETA_PEC)
    phi_miss_ctl = abs(ctl["phi_h"] / ctl["phi_cf_h"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[control] cavity left natural (a void): beta = {ctl['beta']:.6f} "
            f"(void closed form {BETA_VOID:+.2f}), |dbeta| = {delta_beta:.4f} "
            f"vs band {TH15_BETA_BAND:.4%} (gate record 1.4328)"
            f"\n[control] Phi_h = {ctl['phi_h']:.6f} vs Phi_cf,h(PEC) = "
            f"{ctl['phi_cf_h']:.6f}: miss {phi_miss_ctl:.4%} vs band "
            f"{MAT4_BAND:.0%} (gate record 44.2975 %)",
            flush=True,
        )
    assert delta_beta > TH15_BETA_BAND, (
        f"natural-cavity beta {ctl['beta']:.6f} is within the band of 1: not a control"
    )
    assert phi_miss_ctl > MAT4_BAND, (
        f"natural-cavity Phi miss {phi_miss_ctl:.4%} not over {MAT4_BAND:.0%}: "
        "the Phi gate does not resolve the PEC condition"
    )

    # ---- ParaView ------------------------------------------------------------
    e_re, e_mag, sar_ratio = _paraview_fields(msh, anchor["fields"])
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {"E_real": e_re, "E_magnitude": e_mag, "sar_ratio": sar_ratio},
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[paraview] wrote {xdmf_path}"
            f"\n[paraview] colour by `sar_ratio`: the SAR peaks at the poles "
            f"(field axis) and drops below the far field at the equator."
            f"\n\nAll assertions hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
