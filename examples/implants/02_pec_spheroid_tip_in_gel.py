"""Example (`EX-77`): a rod-tip stand-in in the gel -- the slender spheroid's
tip concentrates SAR, and only when the field runs along it.

`IMP-2` step 1 gates a perfectly conducting prolate spheroid (semi-major
c = 6 mm along z, aspect 5, so a = 1.2 mm) cut out of a box of half-extents
(c, c, 2c) as a hole, in `MAT-4`'s gelled-saline host (sigma = 0.57 S/m,
eps_r = 78) at 64 MHz, degree 2, under the imposed quasi-static field E0 zhat
(the closed-form PEC-spheroid trace on the outer walls). The tip field is
E0 / L_z, so the pointwise tip SAR enhancement is (1/L_z)^2; the *tip-ball*
mean Phi_tip = <|E|^2>/E0^2 over the cells within r_t = 2 a^2 / c of either
tip is what the mesh resolves and what the gate asserts. Turning the same
field crosswise (E0 xhat) removes the tip enhancement altogether.

**It asserts, it does not merely render.** Every record, band and helper is
*imported* from the module that gates the capability
(``tests/validation/test_imp2_spheroid_fem.py``), never restated:

* |beta - 1| <= `TH15_BETA_BAND`; |Phi_tip,h / Phi_tip,cf,h - 1| <=
  `MAT4_BAND`; cavity dofs present and max |E| on them < 1e-12 -- each
  printed beside the gate's own W1 record.
* Negative control, asserted on the same cached mesh: the transverse drive
  misses both bands (asserted only as "> band").

Scope: aspect 5, degree 2, quasi-static imposed trace, one c, ``-n 2``. No
rod, no coil, no tissue, no finite k, no heating, no ASTM F2182 claim.

Run it through the example runner (the ``imp:`` group sources the complex
build automatically; a real build raises)::

    ./scripts/run_examples.sh -e imp:2 -n 2 -t 400

Output lands in ``examples/implants/paraview_output``: open
``implants_02_pec_spheroid_tip_in_gel_combined.xdmf`` and colour by
``sar_ratio_axial`` (cut the x-z plane).
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
from fem_em_solver.utils.analytical import prolate_spheroid_depolarisation

# The gated fixture lives in the test that gates `IMP-2` step 1; the runner
# puts only ``src`` on PYTHONPATH, so the repo root goes on sys.path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation.test_imp2_spheroid_fem import (  # noqa: E402
    BOX_HALF_EXTENTS,
    CAVITY_TAG,
    DEGREE,
    E0,
    EPSILON_R_HOST,
    FREQUENCY_HZ,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    SEMI_MAJOR,
    SIGMA_HIGH,
    TH15_BETA_BAND,
    _axial_field,
    _mesh,
    _rss_gib,
    _semi_minor,
    _solve_spheroid,
    _tip_mask,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_02_pec_spheroid_tip_in_gel"
ASPECT = 5.0


def _inv_lz(aspect):
    return 1.0 / prolate_spheroid_depolarisation(np.sqrt(1.0 - 1.0 / aspect**2))


def _dg0_sar(msh, fields):
    """DG0 |E|^2 / E0^2 (= SAR/SAR0 in a homogeneous host) and the CG1 E_real."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar = fem.Function(dg0, name="sar_ratio")
    e = fields.e_complex
    expr = ufl.real(ufl.inner(e, e)) / (E0**2)
    sar.interpolate(fem.Expression(expr, dg0.element.interpolation_points))
    sar.x.scatter_forward()
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_cg = fem.Function(v_cg, name="E_phasor")
    e_cg.interpolate(e)
    e_cg.x.scatter_forward()
    e_re = fem.Function(v_cg, name="E_real")
    e_re.x.array[:] = np.real(e_cg.x.array)
    e_re.x.scatter_forward()
    return sar, e_re


def _dg0_closed_form(msh):
    """DG0 closed-form |E|^2 / E0^2 at the cell midpoints (numpy)."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(dg0, name="sar_ratio_closed_form")
    mids = dg0.tabulate_dof_coordinates()
    e = _axial_field(ASPECT, "continue")(mids)
    f.x.array[:] = np.sum(np.abs(e) ** 2, axis=1) / E0**2
    f.x.scatter_forward()
    return f


def _dg0_copy(msh, src, name):
    dg0 = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(dg0, name=name)
    f.x.array[:] = np.real(src.x.array)
    f.x.scatter_forward()
    return f


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    a = _semi_minor(ASPECT)
    inv_lz = _inv_lz(ASPECT)
    if comm.rank == 0:
        print("=" * 72)
        print("EX-77 -- a rod-tip stand-in in the gel: elongation, and the field along it")
        print("=" * 72)
        print(
            f"\n[geometry] PEC prolate spheroid c = {SEMI_MAJOR * 1e3:.1f} mm, a = "
            f"{a * 1e3:.2f} mm (aspect {ASPECT:.0f}) cut out of a box of half-extents "
            f"{tuple(h * 1e3 for h in BOX_HALF_EXTENTS)} mm as a hole (cavity = facet tag "
            f"{CAVITY_TAG}, outer walls tag {OUTER_BOUNDARY_TAG}), host sigma = {SIGMA_HIGH} "
            f"S/m, eps_r = {EPSILON_R_HOST}, f = {FREQUENCY_HZ / 1e6:.3f} MHz, degree "
            f"{DEGREE}, imposed E0"
            f"\n[closed form] tip field E0/L_z: 1/L_z = {inv_lz:.6f} (gate record "
            f"17.914415, `20261001T201413Z_IMP-2-step1-W0.log:96`); pointwise tip SAR "
            f"enhancement (1/L_z)^2 = {inv_lz**2:.1f} (gate record 320.9)",
            flush=True,
        )
        print("[trend] 1/L_z by aspect (closed form), beside the W0 records:", flush=True)
        for aspect, record in (
            (1.01, "3.024038"),
            (2.0, "5.761564"),
            (5.0, "17.914415"),
            (20.0, "148.168897"),
        ):
            print(
                f"[trend] aspect {aspect:>5}: 1/L_z = {_inv_lz(aspect):.6f} "
                f"(W0 record {record}); (1/L_z)^2 = {_inv_lz(aspect) ** 2:.1f}",
                flush=True,
            )

    # One mesh cache: the gate's own `_mesh()` serves both solves, the figure
    # and the XDMF.
    msh, cell_tags, facet_tags = _mesh(ASPECT)

    # EX-57 setup figure: cell tag 2 only (the spheroid is a hole); clip y < 0
    # opens the block, slice at y = 0 shows the x-z plane with the long axis and
    # the graded tips. No-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={2: "gel host sigma 0.57 S/m (PEC spheroid cut out as a hole)"},
        clip_normal=(0.0, -1.0, 0.0),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:2 -- PEC spheroid (c = 6 mm, aspect 5) cut out of a gel host: x-z slice, tips graded",
        comm=comm,
    )

    # ---- axial anchor (collective: every rank) -------------------------------
    solve_started = time.perf_counter()
    anchor = _solve_spheroid(ASPECT, "axial", return_fields=True)
    sar_axial, e_real_axial = _dg0_sar(msh, anchor.pop("fields"))
    sar_axial_out = _dg0_copy(msh, sar_axial, "sar_ratio_axial")
    sar_cf_out = _dg0_closed_form(msh)
    tip_mask = _tip_mask(msh, ASPECT)
    tip_mask_out = _dg0_copy(msh, tip_mask, "tip_mask")
    e_real_out = e_real_axial
    e_real_out.name = "E_real_axial"
    solve_seconds = time.perf_counter() - solve_started

    beta = anchor["beta"]
    beta_miss = abs(beta - 1.0)
    phi_h, phi_cf_h = anchor["phi_h"], anchor["phi_cf_h"]
    phi_miss = abs(phi_h / phi_cf_h - 1.0)
    rss_axial = _rss_gib(comm)  # collective: outside the rank-0 branch

    if comm.rank == 0:
        print(
            f"\n[solve] axial: {anchor['ncells']} cells (gate record 44 304), tip-mask "
            f"cells {anchor['n_mask']} (6 264), h_tip/rho_t {anchor['h_tip_over_rho']:.3f} "
            f"(0.671), {solve_seconds:.1f} s, RSS(sum) {rss_axial:.2f} GiB"
            f"\n[anchor] (i) fitted beta = {beta:.6f}, |beta - 1| = {beta_miss:.4%} vs "
            f"TH15_BETA_BAND = {TH15_BETA_BAND:.4%} (gate record: beta 0.990979, "
            f"`20261001T201608Z_IMP-2-step1-W1.log:450`); |Im E|/|Re E| = "
            f"{anchor['imag_ratio']:.3e} (6.623e-03)"
            f"\n[anchor] (ii) Phi_tip,h = {phi_h:.6f}, Phi_tip,cf,h = {phi_cf_h:.6f}, "
            f"|Phi_tip,h/Phi_tip,cf,h - 1| = {phi_miss:.4%} vs MAT4_BAND = "
            f"{MAT4_BAND:.0%} (gate record: 49.805748 vs 49.904262, 0.1974 %)"
            f"\n[anchor] pointwise tip (1/L_z)^2 = {inv_lz**2:.1f} vs the tip-ball mean "
            f"{phi_h:.1f}: averaging over r_t = 2 a^2/c = {2 * a * a / SEMI_MAJOR * 1e3:.2f} "
            f"mm takes the peak down by {phi_h / inv_lz**2:.3f}x"
            f"\n[structure] (iii) cavity dofs (reduced) = {anchor['n_cavity']}, max |E| on "
            f"them = {anchor['cavity_max']:.3e} (gate record: 23 886, 0)",
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
        f"Phi_tip,h = {phi_h:.6f} vs Phi_tip,cf,h = {phi_cf_h:.6f}: miss {phi_miss:.4%} "
        f"> {MAT4_BAND:.0%} (imported `MAT-4` band) -- an example/test divergence"
    )

    # ---- negative control: asserted, same mesh, field turned crosswise -------
    ctl = _solve_spheroid(ASPECT, "transverse", return_fields=True)
    sar_tr, _ = _dg0_sar(msh, ctl.pop("fields"))
    sar_tr_out = _dg0_copy(msh, sar_tr, "sar_ratio_transverse")
    delta_beta = abs(ctl["beta"] - 1.0)
    phi_miss_ctl = abs(ctl["phi_h"] / phi_cf_h - 1.0)
    rss_ctl = _rss_gib(comm)
    if comm.rank == 0:
        print(
            f"\n[control] transverse drive (E0 xhat), same mesh ({ctl['ncells']} cells), "
            f"{ctl['seconds']:.1f} s, RSS(sum) {rss_ctl:.2f} GiB: beta = {ctl['beta']:.6f} "
            f"(gate record 0.000044), |beta - 1| = {delta_beta:.4f} vs band "
            f"{TH15_BETA_BAND:.4%}"
            f"\n[control] Phi_tip,h = {ctl['phi_h']:.6f} vs the axial closed form "
            f"{phi_cf_h:.6f}: miss {phi_miss_ctl:.4%} vs band {MAT4_BAND:.0%} (gate "
            f"record 98.5602 %); cavity dofs {ctl['n_cavity']}, max |E| "
            f"{ctl['cavity_max']:.3e}",
            flush=True,
        )
    assert delta_beta > TH15_BETA_BAND, (
        f"transverse beta {ctl['beta']:.6f} is within the band of 1: not a control"
    )
    assert phi_miss_ctl > MAT4_BAND, (
        f"transverse Phi miss {phi_miss_ctl:.4%} not over {MAT4_BAND:.0%}: "
        "the Phi gate does not resolve the orientation"
    )

    # ---- ParaView ------------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {
            "sar_ratio_axial": sar_axial_out,
            "sar_ratio_closed_form": sar_cf_out,
            "sar_ratio_transverse": sar_tr_out,
            "tip_mask": tip_mask_out,
            "E_real_axial": e_real_out,
        },
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[paraview] wrote {xdmf_path}"
            f"\n[paraview] colour by `sar_ratio_axial` (SAR/SAR0, DG0) on an x-z cut: "
            f"the hot spots are the two tips; `sar_ratio_transverse` has none."
            f"\n\nAll assertions hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
