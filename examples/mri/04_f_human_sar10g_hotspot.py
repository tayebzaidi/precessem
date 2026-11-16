"""Example (`EX-62`): the 10 g SAR hotspot at human scale.

`mri:3` / `ans:2` average 10 g SAR on the 4-leg F-small birdcage under a
*single* drive; `ports:15` drives the 16-leg rung in quadrature but reads
only `|B1+|`; `mesh:12` builds the F-human rung without solving it. No
example in the corpus solves a human-scale coil for a safety quantity --
this one does: quadrature drive (ccw/cw superposition) on the refined
F-human phantom, a nine-point candidate lattice for the local-SAR maximum
``c*``, and the C16 / mirror symmetry identities read off the 10 g SAR
itself.

**Import, never restate** (`ANS-1`). Every band, tag and constant --
``C4_COVARIANCE_BAND``, ``EXACT_IDENTITY_RTOL``, ``PHANTOM_CELL_TAG``,
``PHANTOM_RHO_KG_PER_M3``, ``QUADRATURE_DEGREE``, ``TEN_GRAM_KG`` -- comes
from ``tests/validation/test_wf7_f_human_sar10g_hotspot.py``'s own import
chain.

**Rule (a) additive lift, disclosed** (the `EX-61` / `EX-53` precedent). The
gate module's ``sar10g_case`` pytest fixture body is lifted verbatim to a
module-level ``build_sar10g_case(comm)``; the fixture keeps its
``FEM_EM_WF7_STEP3`` env skip and becomes a one-line delegate. No existing
test's behaviour changes -- the fixture still returns exactly what it
returned before, plus additive-only return keys (``mesh``, ``cell_tags``,
``sigma``, ``e_ccw``, ``rho_field``, ``centre``, ``z_c``, ``radius``,
``half_height``, ``a10``, ``az_ref``, ``coords``, ``key_star``, ``c_star``,
``sar_star``, ``control_sar``, ``control``, ``control_worst``,
``phantom_cells``, ``cells``, ``case``, ``build_s``) that an external caller
needs and were not exposed before. The gate module was re-run green from
`main` in this slot (see the guide's harness log for the rerun window).

**Anchors (asserted, all imported):**

* (i) coverage identity -- the whole-phantom ball's dissipated power and
  mass equal the tagged integral, both ``<= EXACT_IDENTITY_RTOL``;
* (ii) the worst of the 15 C16 images of ``c*`` is ``<= C4_COVARIANCE_BAND``;
* (iii) the cw-at-mirror vs ccw-at-``c*`` identity is ``<=`` the same band;
* (iv) every averaging ball (lattice, C16, mirror, control) sits inside the
  measured phantom with 1 mm clear.

**Negative control (asserted -- rule (e), backed by the gate's own
measurement on this fixture).** The single-drive P17 10 g SAR spread over
``c*`` and its three C4 images is asserted ``> C4_COVARIANCE_BAND``; its
ratio to the band and to the ccw worst C16 miss are printed (the gate read
19.28x the band, 285.8x the ccw worst; leg-box-free fixture, GEO-35 2b-ii).

**Printed, never gated:** ``c*`` = r 45 mm, z = z_c; the nine-ball location
table; the whole-phantom mean SAR; step 2's ``|B1+|`` C16/mirror readings on
the same refined mesh; the PRICE line (build / sweep / SAR-phase timing,
summed ``ru_maxrss``).

**Written:** combined XDMF of the pointwise local SAR
``sigma|E|^2/(2*rho)`` on the phantom (ccw drive), plus a DG0 marker field
for the nine lattice balls and a separate marker for ``c*`` (the same
UFL-conditional-ball convention `mri:3` (`EX-53`) uses, applied here to nine
sites instead of four); the nine-ball table as CSV; setup figure per `EX-57`.

**Scope.** One phantom, one fixture, 64 MHz, degree 1. No compliance and no
absolute SAR claim; degree 1 at human scale carries the 5.50% order caveat
(`WF-7` step 0b / 0d).

Needs the complex DolfinX build and ``-n 8`` (heavy tier); the runner sources
the complex build for the ``mri:`` group automatically::

    ./run_examples.sh -e mri:4 -n 8 -t 900

**Outputs** (``paraview_output/``, full filenames -- census references them
this way, never a stem):

* ``mri_04_f_human_sar10g_hotspot_lattice.csv`` -- the nine-ball table (r,
  dz, x, y, z, SAR10g, mass, wall_s, is_c_star);
* ``mri_04_f_human_sar10g_hotspot_combined.xdmf`` (+ ``.h5``) -- pointwise
  ``SAR`` (DG0), ``LatticeBall`` (DG0, 1-9 = the nine lattice sites),
  ``CStar`` (DG0, 1 at the argmax site) and ``CellTags``; plus the ccw drive's
  fields over the whole domain (added 2026-09-30, operator request):
  ``E_real`` / ``E_imag`` / ``E_magnitude`` (CG1), ``B_real`` / ``B_imag`` /
  ``B_magnitude`` and ``B1_plus`` / ``B1_minus`` (DG0) — printed only.

Setup figure per `EX-57` in ``figures/mri_04_f_human_sar10g_hotspot_setup.png``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.faraday import (  # noqa: E402
    b1_minus,
    b1_plus,
    magnetic_flux_density_from_e,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_birdcage_b1_plus_map import (  # noqa: E402
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ  # noqa: E402
from tests.validation.test_wf7_f_human_sar10g_hotspot import (  # noqa: E402
    build_sar10g_case,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_04_f_human_sar10g_hotspot"


def _pointwise_sar_field(mesh, e_complex, sigma_field):
    """``sigma(x)|E(x)|^2 / (2 rho)`` as a DG0 field, ``rho`` the phantom's
    physical density -- `mri:3`'s (`EX-53`) route: a **scalar**
    ``PHANTOM_RHO_KG_PER_M3``, never the gate's zero-elsewhere ``rho_field``
    (that field exists only so a ball's *mass* integral excludes
    non-phantom cells; used pointwise it would divide by zero everywhere
    outside the phantom)."""
    expr = 0.5 * sigma_field * ufl.inner(e_complex, e_complex) / PHANTOM_RHO_KG_PER_M3
    dg0 = fem.functionspace(mesh, ("DG", 0))
    sar = fem.Function(dg0, name="SAR")
    sar.interpolate(fem.Expression(expr, dg0.element.interpolation_points))
    sar.x.array[:] = np.real(sar.x.array)
    sar.x.scatter_forward()
    return sar


def _ball_marker_field(mesh, centres, radius, name):
    """DG0 integer field: 0 outside every ball, ``k + 1`` inside ball ``k`` --
    the same UFL ``conditional`` on squared distance that ``mass_averaged_sar``
    uses internally (`post/sar.py`), the `mri:3` convention applied to
    ``len(centres)`` sites instead of four."""
    x = ufl.SpatialCoordinate(mesh)
    value = ufl.as_ufl(0.0)
    for k, c in enumerate(centres):
        offset = x - ufl.as_vector([float(v) for v in c])
        inside = ufl.lt(ufl.real(ufl.dot(offset, offset)), float(radius) ** 2)
        value = ufl.conditional(inside, float(k + 1), value)

    dg0 = fem.functionspace(mesh, ("DG", 0))
    field = fem.Function(dg0, name=name)
    field.interpolate(fem.Expression(value, dg0.element.interpolation_points))
    field.x.array[:] = np.real(field.x.array)
    field.x.scatter_forward()
    return field


def _field_export(e_complex, omega):
    """The ccw drive's E and B for ParaView — whole domain, not only the phantom.

    XDMF carries no complex arrays, so each phasor is split into real / imag
    vector fields (the `ports:4` / `ports:7` convention): ``E`` interpolated
    from N1curl to CG1 vectors, ``B = curl E / (-jω)`` as DG0 vectors
    (:func:`magnetic_flux_density_from_e`), their magnitudes, and the rotating
    components ``|B1+|`` / ``|B1-|`` (DG0, :func:`b1_plus` / :func:`b1_minus`,
    real magnitudes stored in complex arrays, hence ``.real``). Peak-phasor
    convention, same drive as ``SAR``. Printed only — nothing here is gated.
    """
    msh = e_complex.function_space.mesh
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_cg = fem.Function(v_cg)
    e_cg.interpolate(e_complex)
    e_cg.x.scatter_forward()
    b_dg = magnetic_flux_density_from_e(e_complex, omega)

    s_cg = fem.functionspace(msh, ("Lagrange", 1))
    s_dg = fem.functionspace(msh, ("DG", 0))
    out = {}
    for label, src, scalar_space in (("E", e_cg, s_cg), ("B", b_dg, s_dg)):
        re = fem.Function(src.function_space, name=f"{label}_real")
        re.x.array[:] = np.real(src.x.array)
        im = fem.Function(src.function_space, name=f"{label}_imag")
        im.x.array[:] = np.imag(src.x.array)
        mag = fem.Function(scalar_space, name=f"{label}_magnitude")
        comps = np.abs(src.x.array.reshape(-1, 3))
        mag.x.array[:] = np.sqrt(np.sum(comps * comps, axis=1))
        out.update({f"{label}_real": re, f"{label}_imag": im, f"{label}_magnitude": mag})
    for name, fn in (("B1_plus", b1_plus(b_dg)), ("B1_minus", b1_minus(b_dg))):
        real = fem.Function(fn.function_space, name=name)
        real.x.array[:] = np.real(fn.x.array)
        out[name] = real
    for f in out.values():
        f.x.scatter_forward()
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mri:` "
            "group sources it automatically"
        )

    if comm.rank == 0:
        print("=" * 78, flush=True)
        print(
            "EX-62 -- the 10 g SAR hotspot at human scale (WF-7 step 3, "
            "quadrature drive)",
            flush=True,
        )
        print("=" * 78, flush=True)
        print(
            "\n[scope] one phantom, one fixture, 64 MHz, degree 1. No compliance, "
            "no absolute SAR claim; degree 1 at human scale carries the 5.50% "
            "order caveat (WF-7 step 0b/0d).",
            flush=True,
        )
        print(
            f"\n[gate import] build_sar10g_case from "
            f"tests/validation/test_wf7_f_human_sar10g_hotspot.py "
            f"(rule (a) additive lift, disclosed above); -n {comm.size}",
            flush=True,
        )

    # ---- run the gate's own construction, verbatim -------------------------
    case = build_sar10g_case(comm)

    if comm.rank == 0:
        print(
            f"\n[gate reproduced] {case['cells']} cells "
            f"({case['phantom_cells']} phantom cells), build {case['build_s']:.2f} s, "
            f"window {case['t_window']:.2f} s at -n {comm.size}",
            flush=True,
        )

    # EX-57 setup figure -- right after the fixture returns the built mesh.
    write_setup_figure(
        case["mesh"],
        case["cell_tags"],
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={CONDUCTOR_CELL_TAG: "conductor", PHANTOM_CELL_TAG: "phantom (saline)"},
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, float(case["z_c"])),
        title="mri:4 -- F-human 10 g SAR hotspot (WF-7 step 3, ccw drive)",
        comm=comm,
    )

    # ---- (iv) containment, printed then asserted ----------------------------
    if comm.rank == 0:
        print(
            f"\n[containment] {len(case['containment'])} balls checked, 1 mm "
            "clear (ASSERTED):",
            flush=True,
        )
        worst_excess = max(e for _lab, _c, e in case["containment"])
        print(f"    worst excess {worst_excess:.6e} m (<= 0.0 means contained)", flush=True)
    bad = [(lab, excess) for lab, _c, excess in case["containment"] if excess > 0.0]
    assert not bad, f"balls outside the phantom (1 mm margin): {bad}"

    # ---- (i) coverage identity ----------------------------------------------
    if comm.rank == 0:
        print(
            f"\n[anchor i] coverage ball radius {case['coverage_radius'] * 1e3:.4f} mm "
            f"(clear of the conductor: < {case['conductor_limit'] * 1e3:.1f} mm): "
            f"power rel {case['power_miss']:.3e}, mass rel {case['mass_miss']:.3e} "
            f"(ASSERTED both <= EXACT_IDENTITY_RTOL {case['rtol']:.0e})",
            flush=True,
        )
    assert case["coverage_radius"] < case["conductor_limit"], (
        f"coverage ball {case['coverage_radius']:.6f} m reaches within the "
        "conductor clearance"
    )
    assert case["power_miss"] <= case["rtol"], f"power {case['power_miss']:.3e} > {case['rtol']:g}"
    assert case["mass_miss"] <= case["rtol"], f"mass {case['mass_miss']:.3e} > {case['rtol']:g}"

    # ---- (ii) C16 -------------------------------------------------------------
    c16 = case["c16"]
    worst_m = case["worst_m"]
    if comm.rank == 0:
        print(
            f"\n[anchor ii] worst of the 15 C16 images: {c16[worst_m] * 100:.4f}% "
            f"at R_{worst_m} (ASSERTED <= C4_COVARIANCE_BAND "
            f"{case['band'] * 100:.1f}%)",
            flush=True,
        )
    assert c16[worst_m] <= case["band"], (
        f"10 g SAR not C16-invariant: R_{worst_m} {c16[worst_m] * 100:.4f}% > "
        f"{case['band'] * 100:.1f}%"
    )

    # ---- (iii) mirror -----------------------------------------------------------
    if comm.rank == 0:
        print(
            f"\n[anchor iii] mirror {case['mirror'] * 100:.4f}% (ASSERTED <= "
            f"{case['band'] * 100:.1f}%)",
            flush=True,
        )
    assert case["mirror"] <= case["band"], (
        f"mirror {case['mirror'] * 100:.4f}% > {case['band'] * 100:.1f}%"
    )

    # ---- negative control: single-drive P17 spread over C4 images -----------
    control_worst = case["control_worst"]
    band_ratio = control_worst / case["band"]
    ccw_ratio = control_worst / max(c16[worst_m], 1e-300)
    if comm.rank == 0:
        print(
            f"\n[negative control] single-drive P17 10 g SAR spread over c* and "
            f"its three C4 images: {control_worst * 100:.4f}% "
            f"(ASSERTED > C4_COVARIANCE_BAND {case['band'] * 100:.1f}%); "
            f"{band_ratio:.4f}x the band, {ccw_ratio:.1f}x the ccw worst "
            f"({c16[worst_m] * 100:.4f}%)",
            flush=True,
        )
    assert control_worst > case["band"], (
        f"single-drive P17 spread {control_worst * 100:.4f}% is inside the "
        f"{case['band'] * 100:.1f}% band -- the quadrature drive is not doing "
        "the work this control is supposed to isolate"
    )

    # ---- printed only: c*, the nine-ball table, step 2's |B1+| readings -----
    reads = case["reads"]
    coords = case["coords"]
    key_star = case["key_star"]
    if comm.rank == 0:
        print(
            f"\n[c*] r {key_star[0] * 1e3:.0f} mm, z = z_c{key_star[1] * 1e3:+.0f} mm "
            f"({case['c_star'][0][0] * 1e3:.3f}, {case['c_star'][0][1] * 1e3:.3f}, "
            f"{case['c_star'][0][2] * 1e3:.3f}) mm, SAR10g(c*) "
            f"{case['sar_star']:.9e} W/kg",
            flush=True,
        )
        print("\n[nine-ball location table, printed]", flush=True)
        for (r, dz), out in sorted(reads.items(), key=lambda kv: (kv[0][1], kv[0][0])):
            marker = "   <- c*" if (r, dz) == key_star else ""
            print(
                f"    r = {r * 1e3:4.0f} mm, z = z_c{dz * 1e3:+5.0f} mm: SAR10g "
                f"{out['averaged_sar_w_per_kg']:.9e} W/kg, mass "
                f"{out['mass_kg'] * 1e3:.6f} g, {out['wall_s']:.2f} s{marker}",
                flush=True,
            )
        b1 = case["case"]
        print(
            f"\n[step 2 |B1+|, printed, same refined mesh] worst C16 "
            f"{b1['c16'][b1['worst_m']] * 100:.4f}% at R_{b1['worst_m']}, mirror "
            f"{b1['mirror'] * 100:.4f}%, {b1['n_valid']}/{b1['points'].shape[0]} "
            f"valid sample points, exact identity {b1['identity_dev']:.3e}",
            flush=True,
        )
        print(
            f"\n[PRICE] build {case['build_s']:.2f} s, window "
            f"{case['t_window']:.2f} s wall at -n {comm.size}; summed ru_maxrss "
            f"{case['rss_after_sweep']:.3f} GiB after the sweep, "
            f"{case['rss_final']:.3f} GiB at the end",
            flush=True,
        )

    # ---- CSV of the nine-ball table -------------------------------------------
    if comm.rank == 0:
        OUTPUT_DIR.mkdir(exist_ok=True)
        csv_path = OUTPUT_DIR / f"{BASENAME}_lattice.csv"
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                ["r_m", "dz_m", "x_m", "y_m", "z_m", "sar10g_w_per_kg", "mass_kg", "wall_s", "is_c_star"]
            )
            for (r, dz), out in sorted(reads.items(), key=lambda kv: (kv[0][1], kv[0][0])):
                c = coords[(r, dz)]
                writer.writerow(
                    [
                        f"{r:.6f}",
                        f"{dz:.6f}",
                        f"{c[0]:.9e}",
                        f"{c[1]:.9e}",
                        f"{c[2]:.9e}",
                        f"{out['averaged_sar_w_per_kg']:.9e}",
                        f"{out['mass_kg']:.9e}",
                        f"{out['wall_s']:.3f}",
                        "1" if (r, dz) == key_star else "0",
                    ]
                )
        adopt_host_ownership(OUTPUT_DIR, comm=comm)
        print(f"\n[csv] {csv_path}", flush=True)

    # ---- ParaView: pointwise SAR + nine-ball / c* markers ----------------------
    lattice_order = sorted(reads.keys(), key=lambda k: (k[1], k[0]))
    lattice_centres = [coords[k] for k in lattice_order]
    sar_field = _pointwise_sar_field(case["mesh"], case["e_ccw"], case["sigma"])
    lattice_field = _ball_marker_field(
        case["mesh"], lattice_centres, case["a10"], "LatticeBall"
    )
    c_star_field = _ball_marker_field(
        case["mesh"], [tuple(float(v) for v in case["c_star"][0])], case["a10"], "CStar"
    )

    field_export = _field_export(case["e_ccw"], 2.0 * np.pi * FREQUENCY_HZ)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        case["mesh"],
        case["cell_tags"],
        {"SAR": sar_field, "LatticeBall": lattice_field, "CStar": c_star_field,
         **field_export},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[paraview] wrote {xdmf_path}"
            "\n[paraview] threshold `LatticeBall` (1-9) to see the nine candidate "
            "sites; `CStar` (1) marks the argmax alone; colour by `SAR` (W/kg)."
            "\n[paraview] fields of the same ccw drive, whole domain: `E_real` / "
            "`E_imag` / `E_magnitude` (CG1, V/m), `B_real` / `B_imag` / "
            "`B_magnitude` (DG0, T), `B1_plus` / `B1_minus` (DG0, T) — peak phasors."
            f"\n\nAll asserted anchors hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
