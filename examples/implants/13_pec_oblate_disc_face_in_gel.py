"""Example (`EX-147`): a PEC disc face-on to the field -- the pole enhancement
1/L_z of an oblate spheroid in gel, and the edge-on drive that misses it.

`IMP-6` step 3 gates the **axial** (face-on) drive of a perfectly conducting
oblate spheroid (a = 6 mm in x and y, c = a/5 = 1.2 mm along z, aspect 5) cut
out of a gelled-saline host (sigma = 0.57 S/m, eps_r = 78) at 64 MHz, degree 2,
under the imposed quasi-static field E0 zhat.  The field at the flat faces is
E0 / L_z, so the pole-ball mean Phi_pole = <|E|^2> / E0^2 over the cells within
r_p = a/5 of (0, 0, +-c) is what the mesh resolves and what the gate asserts.
`imp:8` shows the same disc edge-on (the rim, 1/L_x); this example is its
face-on twin.

**It asserts, it does not merely render.**  Every record, band and helper is
*imported* from the modules that gate the capability, never restated:

* the least-squares amplitude beta_axial: |beta - 1| <= `TH15_BETA_BAND`;
* |Phi_pole,h / Phi_pole,cf,h - 1| <= `MAT4_BAND`;
* cavity dofs present and max |E| on them below the gate's own literal
  (all three inside `_report_and_gate`, the gate's own function);
* negative control on the same cached mesh: the in-plane drive (E0 xhat)
  misses both bands (asserted only as "> band").

Scope: aspect 5, a = 6 mm, degree 2, quasi-static imposed trace, one mesh,
``-n 2``.  Aspect 2 (gate W1a) is not run.  No coil, no tissue, no SAR
averaging, no heating, no external check.

Run it through the example runner (the ``imp:`` group sources the complex
build automatically; a real build raises)::

    ./run_examples.sh -e imp:13 -n 2 -t 600

Output lands in ``examples/implants/paraview_output``: open
``implants_13_pec_oblate_disc_face_in_gel_combined.xdmf`` and colour by
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
from fem_em_solver.utils.analytical import oblate_spheroid_depolarisation

# The gated fixtures live in the tests that gate `IMP-6` step 3 / step 1; the
# runner puts only ``src`` on PYTHONPATH, so the repo root goes on sys.path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_imp6_oblate_axial_fem as m6a  # noqa: E402
import tests.validation.test_imp6_oblate_fem as m6  # noqa: E402
from tests.validation.test_imp1_pec_sphere_lossy_host import (  # noqa: E402
    E0,
    MAT4_BAND,
    TH15_BETA_BAND,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_13_pec_oblate_disc_face_in_gel"
ASPECT = 5.0

# Records the example prints beside its readings (printed, never asserted by
# equality).  Log:line of each is cited in the guide.
REC = dict(
    beta_axial=1.005469,  # 20261013T050558Z_IMP-6-step3-W1b.log:190
    imag_ratio=9.846e-03,  # :190
    phi_h=1.683758,  # :191
    phi_cf_h=1.667975,  # :191
    miss_pct=0.9463,  # :191
    n_mask=154,  # :189
    h_pole_mm=1.948,  # :189
    deficit=-5.777e-03,  # :189
    n_cavity=19952,  # :192
    ctl_beta=-0.032226,  # 20261013T050819Z_IMP-6-step3-W2.log:192
    ctl_phi_h=0.043756,  # :193
    ctl_miss_pct=97.3767,  # :193
    l_z=0.750484,  # W0 :303
)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    _, l_z = oblate_spheroid_depolarisation(float(np.sqrt(1.0 - 1.0 / ASPECT**2)))
    inv_lz = 1.0 / l_z

    if comm.rank == 0:
        print("=" * 72)
        print("EX-147 -- a PEC disc face-on to the field: the pole, and the drive that misses it")
        print("=" * 72)
        print(
            f"\n[closed form] L_z = {l_z:.6f} (record {REC['l_z']}); flat-face field "
            f"E0/L_z: 1/L_z = {inv_lz:.6f}; pointwise (1/L_z)^2 = {inv_lz**2:.3f}; pole ball "
            f"r_p = a/5 = {m6a.POLE_RADIUS * 1e3:.2f} mm",
            flush=True,
        )

    # ---- gate solve (collective: every rank) ----------------------------------
    d = m6a._solve_pole(ASPECT, "axial")
    if comm.rank == 0:
        print(
            f"\n[records] W1b: beta_axial {REC['beta_axial']}, |Im E|/|Re E| {REC['imag_ratio']:.3e}, "
            f"Phi_pole,h {REC['phi_h']} vs Phi_pole,cf,h {REC['phi_cf_h']} (+{REC['miss_pct']} %), "
            f"pole-mask cells {REC['n_mask']}, h_pole {REC['h_pole_mm']} mm, deficit "
            f"{REC['deficit']:.3e}, cavity dofs {REC['n_cavity']} "
            f"(`20261013T050558Z_IMP-6-step3-W1b.log:189-192`); gate prints and asserts follow",
            flush=True,
        )
    m6a._report_and_gate(ASPECT, d)  # beta, Phi, cavity asserts (imported)
    phi_h, phi_cf_h = d["phi_h"], d["phi_cf_h"]
    if comm.rank == 0:
        print(
            f"[anchor] ball mean Phi_pole,cf,h {phi_cf_h:.6f} sits below the pointwise "
            f"(1/L_z)^2 = {inv_lz**2:.3f}: the ball reaches off the pole "
            f"(ratio {phi_cf_h / inv_lz**2:.3f})",
            flush=True,
        )

    # ---- negative control: in-plane drive, same mesh --------------------------
    c = m6a._solve_pole(ASPECT, "inplane")
    c_beta_miss = abs(c["beta"] - 1.0)
    c_phi_miss = abs(c["phi_h"] / c["phi_cf_h"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[control] in-plane drive (E0 xhat), same mesh ({c['ncells']} cells), "
            f"{c['seconds']:.1f} s: beta (axial basis) = {c['beta']:.6f} (record "
            f"{REC['ctl_beta']}), |beta - 1| = {c_beta_miss:.4f} vs band {TH15_BETA_BAND:.4%} "
            f"({c_beta_miss / TH15_BETA_BAND:.1f}x)"
            f"\n[control] Phi_pole,h = {c['phi_h']:.6f} (record {REC['ctl_phi_h']}) vs axial "
            f"Phi_pole,cf,h {c['phi_cf_h']:.6f}: miss {c_phi_miss:.4%} vs band "
            f"{MAT4_BAND:.0%} (record {REC['ctl_miss_pct']} %, {c_phi_miss / MAT4_BAND:.1f}x)",
            flush=True,
        )
    assert c_beta_miss > TH15_BETA_BAND, (
        f"in-plane beta {c['beta']:.6f} is within the band of 1: the drive did not switch"
    )
    assert c_phi_miss > MAT4_BAND, (
        f"in-plane Phi miss {c_phi_miss:.4%} not over {MAT4_BAND:.0%}: the drive did not switch"
    )

    # ---- picture solve (collective) -------------------------------------------
    msh, cell_tags, facet_tags = m6._mesh(ASPECT)
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={2: "gel host sigma 0.57 S/m (PEC oblate spheroid cut out as a hole)"},
        clip_normal=(0.0, -1.0, 0.0),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:13 -- PEC oblate disc (a = 6 mm, aspect 5) face-on to E0: x-z slice, poles at z = +-c",
        comm=comm,
    )
    pic = m6._solve_oblate(ASPECT, "axial", return_fields=True)
    fields = pic.pop("fields")
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar_ax = fem.Function(dg0, name="sar_ratio_axial")
    e = fields.e_complex
    sar_ax.interpolate(
        fem.Expression(ufl.real(ufl.inner(e, e)) / (E0**2), dg0.element.interpolation_points)
    )
    sar_ax.x.scatter_forward()
    pole_mask = m6a._pole_mask(msh, ASPECT)
    pole_mask.name = "pole_mask"
    del fields, e

    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {"sar_ratio_axial": sar_ax, "pole_mask": pole_mask},
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[picture] same-mesh axial solve: {pic['ncells']} cells, {pic['seconds']:.1f} s; "
            f"beta (in-plane basis, the routine's fit; record of imp:8's control) {pic['beta']:.6f}"
            f"\n[paraview] wrote {xdmf_path}"
            f"\n[paraview] colour by `sar_ratio_axial` (|E|^2/E0^2, DG0) on an x-z cut: the "
            f"faces at z = +-c carry the pole enhancement, the rim stays cold; `pole_mask` "
            f"marks the pole-ball cells."
            f"\n[orientation] face-on, a thin disc concentrates the field at its faces only by "
            f"1/L_z = {inv_lz:.3f}; edge-on (`imp:8`) its rim reaches 1/L_x = 8.016."
            f"\n\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
