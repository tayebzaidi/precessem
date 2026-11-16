"""Example (`EX-111`): a disc-edge stand-in in the gel -- a flat PEC body's rim
concentrates the field only under an in-plane drive, and less than a rod tip of
the same curvature.

`IMP-6` step 1 gates a perfectly conducting **oblate** spheroid (equatorial
semi-axis a = 6 mm in x and y, polar c = a/5 = 1.2 mm along z, aspect 5) cut out
of a box of half-extents (1.75a, 1.75a, a) as a hole, in `MAT-4`'s gelled-saline
host (sigma = 0.57 S/m, eps_r = 78) at 64 MHz, degree 2, under the imposed
quasi-static **in-plane** field E0 xhat (the closed-form PEC-oblate-spheroid
trace on the outer walls).  The rim field is E0 / L_x, so the pointwise rim SAR
enhancement is (1/L_x)^2; the *rim-ball* mean Phi_rim = <|E|^2>/E0^2 over the
cells within r_r = 2 c^2 / a of either rim point (+-a, 0, 0) is what the mesh
resolves and what the gate asserts.  Driven face-on (E0 zhat) the same disc
leaves its rim cold.

**It asserts, it does not merely render.**  Every record, band and helper is
*imported* from the module that gates the capability
(``tests/validation/test_imp6_oblate_fem.py``), never restated:

* |beta - 1| <= `TH15_BETA_BAND`; |Phi_rim,h / Phi_rim,cf,h - 1| <= `MAT4_BAND`;
  cavity dofs present and max |E| on them < 1e-12 -- each printed beside the
  gate's own W1b record.
* Negative control, asserted on the same cached mesh: the axial (face-on)
  drive misses both bands (asserted only as "> band").

Scope: aspect 5, a = 6 mm, degree 2, quasi-static imposed trace, one mesh,
``-n 2``.  A smooth spheroid, not a square-edged plate (a real plate edge
carries a wedge singularity the spheroid does not).  No coil, no tissue, no
finite k, no SAR averaging, no heating, no ASTM F2182 claim.

Run it through the example runner (the ``imp:`` group sources the complex
build automatically; a real build raises)::

    ./scripts/run_examples.sh -e imp:8 -n 2 -t 400

Output lands in ``examples/implants/paraview_output``: open
``implants_08_pec_oblate_disc_edge_in_gel_combined.xdmf`` and colour by
``sar_ratio_inplane`` (cut the x-z plane); the rim-line CSV is
``implants_08_pec_oblate_disc_edge_in_gel_rim_line.csv``.
"""

from __future__ import annotations

import math
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
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.utils.analytical import prolate_spheroid_depolarisation

# The gated fixtures live in the tests that gate `IMP-6` step 1 / `IMP-2`
# step 1; the runner puts only ``src`` on PYTHONPATH, so the repo root goes on
# sys.path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_imp2_spheroid_fem as m2  # noqa: E402
import tests.validation.test_imp6_oblate_fem as m6  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_08_pec_oblate_disc_edge_in_gel"
ASPECT = 5.0
E0 = m6.E0

# Records the example prints beside its readings (printed, never asserted by
# equality).  Every entry carries the log:line it was read from.
REC = dict(
    ncells=32676,  # 20261012T184537Z_IMP-6-step1-W1b.log:189
    n_mask=4905,  # :189
    h_rim_over_rho=0.665,  # :189
    deficit=-5.777e-03,  # :189
    beta=1.006343,  # :190
    imag_ratio=8.516e-03,  # :190
    phi_h=23.024938,  # :191
    phi_cf_h=22.702995,  # :191
    n_cavity=19952,  # :192
    min_sicn=0.3707,  # :187
    ctl_beta=-0.001796,  # 20261012T184748Z_IMP-6-step1-W2.log:190
    ctl_beta_axial=1.005469,  # :190
    ctl_phi_h=0.491473,  # :191
    ctl_miss=97.8352,  # percent, :191
    inv_lx=8.015515,  # 20261012T183435Z_IMP-6-step0-W1.log:53
    inv_lz_obl=1.332474,  # :53
    inv_lz_pro=17.914415,  # 20261001T201413Z_IMP-2-step1-W0.log:96
    phi_tip_h=49.805748,  # IMP-2 step 1 W1 (EX-77 guide)
)
TREND_LX = {1.01: "3.0120", 2.0: "4.2301", 5.0: "8.0155", 20.0: "27.0935", 100.0: "128.9466"}
TREND_LZ = {1.01: "3.024038", 2.0: "5.761564", 5.0: "17.914415", 20.0: "148.168897"}


def _ecc(aspect):
    return math.sqrt(1.0 - 1.0 / aspect**2)


def _inv_lx(aspect):
    lx, _ = m6.oblate_spheroid_depolarisation(_ecc(aspect))
    return 1.0 / lx


def _inv_lz_pro(aspect):
    return 1.0 / prolate_spheroid_depolarisation(_ecc(aspect))


def _dg0_from(msh, values, name):
    dg0 = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(dg0, name=name)
    f.x.array[:] = np.real(values)
    f.x.scatter_forward()
    return f


def _dg0_sar(msh, fields, name):
    """DG0 |E|^2 / E0^2 (= SAR/SAR0 in a homogeneous host) and a CG1 E_real."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar = fem.Function(dg0, name=name)
    e = fields.e_complex
    expr = ufl.real(ufl.inner(e, e)) / (E0**2)
    sar.interpolate(fem.Expression(expr, dg0.element.interpolation_points))
    sar.x.scatter_forward()
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_cg = fem.Function(v_cg)
    e_cg.interpolate(e)
    e_cg.x.scatter_forward()
    e_re = fem.Function(v_cg, name="E_real_inplane")
    e_re.x.array[:] = np.real(e_cg.x.array)
    e_re.x.scatter_forward()
    return sar, e_re


def _evaluate_abs_e(fields, points, comm):
    """FEM |E| and Re E_x at ``points`` (collective); raises on a missed point."""
    re_vals, ok_re = evaluate_vector_field_parallel(fields.e_real, points, comm)
    im_vals, ok_im = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    ok = np.asarray(ok_re & ok_im)
    if not np.all(ok):
        bad = points[~ok][0]
        raise RuntimeError(f"point {bad} was not evaluated")
    re3 = np.real(re_vals[:, :3])
    im3 = np.real(im_vals[:, :3])
    return np.sqrt(np.sum(re3**2 + im3**2, axis=1)), re3[:, 0]


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    # Structural asserts: the two gate modules share their bands and the 6 mm
    # long semi-axis, and the matched curvature radius the comparison rests on.
    assert m6.TH15_BETA_BAND is m2.TH15_BETA_BAND
    assert m6.MAT4_BAND is m2.MAT4_BAND
    assert m6.SEMI_MAJOR == m2.SEMI_MAJOR == 0.006
    assert math.isclose(
        m6._rho_r(5.0), m2._semi_minor(5.0) ** 2 / m2.SEMI_MAJOR, rel_tol=1e-12
    )

    a = m6.SEMI_MAJOR
    c = m6._semi_minor(ASPECT)
    rho_r = m6._rho_r(ASPECT)
    c_pro = m2.SEMI_MAJOR
    a_pro = m2._semi_minor(ASPECT)
    rho_t = a_pro**2 / c_pro
    inv_lx = _inv_lx(ASPECT)
    lx = 1.0 / inv_lx
    inv_lz_obl = 1.0 / m6.oblate_spheroid_depolarisation(_ecc(ASPECT))[1]
    inv_lz_pro = _inv_lz_pro(ASPECT)

    if comm.rank == 0:
        print("=" * 72)
        print("EX-111 -- a disc-edge stand-in in the gel: the rim, and the field across it")
        print("=" * 72)
        print(
            f"\n[geometry] PEC oblate spheroid a = {a * 1e3:.1f} mm (x, y), c = {c * 1e3:.2f} mm "
            f"(z), aspect {ASPECT:.0f}, cut out of a box of half-extents "
            f"{tuple(h * 1e3 for h in m6.BOX_HALF_EXTENTS)} mm as a hole (cavity = facet tag "
            f"{m6.CAVITY_TAG}, outer walls tag {m6.OUTER_BOUNDARY_TAG}), host sigma = "
            f"{m6.SIGMA_HIGH} S/m, eps_r = {m6.EPSILON_R_HOST}, f = "
            f"{m6.FREQUENCY_HZ / 1e6:.3f} MHz, degree {m6.DEGREE}, imposed E0"
            f"\n[closed form] rim field E0/L_x: L_x = {lx:.6f}, 1/L_x = {inv_lx:.6f} (record "
            f"{REC['inv_lx']}, `20261012T183435Z_IMP-6-step0-W1.log:53`); pointwise rim "
            f"(1/L_x)^2 = {inv_lx**2:.2f}; flat-face field 1/L_z = {inv_lz_obl:.6f} (record "
            f"{REC['inv_lz_obl']})",
            flush=True,
        )
        print("[trend] 1/L_x (oblate rim) and 1/L_z (prolate tip) by aspect, closed form:", flush=True)
        for asp in (1.01, 2.0, 5.0, 20.0, 100.0):
            rec_x = TREND_LX[asp]
            rec_z = TREND_LZ.get(asp, "no record")
            print(
                f"[trend] aspect {asp:>6}: 1/L_x^obl = {_inv_lx(asp):.4f} (step-0 record "
                f"{rec_x}); 1/L_z^pro = {_inv_lz_pro(asp):.6f} (W0 record {rec_z})",
                flush=True,
            )

    # One mesh cache: the gate's own `_mesh()` serves both solves, the figure
    # and the XDMF.
    msh, cell_tags, facet_tags = m6._mesh(ASPECT)

    # EX-57 setup figure: cell tag 2 only (the spheroid is a hole); clip y < 0
    # opens the block, slice at y = 0 shows the x-z plane with the disc edge-on
    # and both rims graded.  No-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={2: "gel host sigma 0.57 S/m (PEC oblate spheroid cut out as a hole)"},
        clip_normal=(0.0, -1.0, 0.0),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:8 -- PEC oblate spheroid (a = 6 mm, aspect 5) cut out of a gel host: x-z slice, rim graded",
        comm=comm,
    )

    # ---- in-plane anchor (collective: every rank) -----------------------------
    solve_started = time.perf_counter()
    d = m6._solve_oblate(ASPECT, "inplane", return_fields=True)
    fields = d.pop("fields")
    sar_in, e_real_out = _dg0_sar(msh, fields, "sar_ratio_inplane")

    dg0 = fem.functionspace(msh, ("DG", 0))
    mids = dg0.tabulate_dof_coordinates()
    e_cf = m6._field(ASPECT, "inplane", "continue")(mids)
    sar_cf_out = _dg0_from(msh, np.sum(np.abs(e_cf) ** 2, axis=1) / E0**2, "sar_ratio_closed_form")
    sar_diff_out = _dg0_from(
        msh, np.real(sar_in.x.array) - np.real(sar_cf_out.x.array), "sar_ratio_minus_closed_form"
    )
    rim_mask_out = _dg0_from(msh, m6._rim_mask(msh, ASPECT).x.array, "rim_mask")

    # ---- samples read from the in-plane field before it is freed --------------
    # Four azimuths at r = a + rho_r/2, z = 0.
    r_az = a + 0.5 * rho_r
    phis = np.deg2rad([0.0, 30.0, 60.0, 90.0])
    az_pts = np.column_stack([r_az * np.cos(phis), r_az * np.sin(phis), np.zeros(4)])
    az_fem, _ = _evaluate_abs_e(fields, az_pts, comm)
    az_cf = np.linalg.norm(np.real(m6._field(ASPECT, "inplane")(az_pts)), axis=1)

    # Rim line: 200 points on y = z = 0, x = a + s, s log-spaced in [rho_r/10, 0.7a].
    s = np.logspace(math.log10(rho_r / 10.0), math.log10(0.7 * a), 200)
    line_pts = np.column_stack([a + s, np.zeros_like(s), np.zeros_like(s)])
    line_fem, line_ex_re = _evaluate_abs_e(fields, line_pts, comm)
    line_cf_vec = np.real(m6._field(ASPECT, "inplane")(line_pts))
    line_cf = np.linalg.norm(line_cf_vec, axis=1)
    pro_pts = np.column_stack([np.zeros_like(s), np.zeros_like(s), c_pro + s])
    line_pro = np.linalg.norm(np.real(m2._axial_field(ASPECT)(pro_pts)), axis=1)
    rel = line_fem / line_cf - 1.0

    del fields  # cap RSS before the control solve
    solve_seconds = time.perf_counter() - solve_started

    beta = d["beta"]
    beta_miss = abs(beta - 1.0)
    phi_h, phi_cf_h = d["phi_h"], d["phi_cf_h"]
    phi_miss = abs(phi_h / phi_cf_h - 1.0)
    rss_in = m6._rss_gib(comm)  # collective: outside the rank-0 branch

    if comm.rank == 0:
        print(
            f"\n[solve] in-plane: {d['ncells']} cells (W1b record {REC['ncells']}), rim-mask cells "
            f"{d['n_mask']} ({REC['n_mask']}), h_rim/rho_r {d['h_rim_over_rho']:.3f} "
            f"({REC['h_rim_over_rho']}), body-volume deficit {d['deficit']:.3e} "
            f"({REC['deficit']:.3e}), {solve_seconds:.1f} s, RSS(sum) {rss_in:.2f} GiB"
            f"\n[record] mesher min SICN printed above (W1b record {REC['min_sicn']}, `:187`)"
            f"\n[anchor] (i) fitted beta = {beta:.6f}, |beta - 1| = {beta_miss:.4%} vs "
            f"TH15_BETA_BAND = {m6.TH15_BETA_BAND:.4%} (gate record: beta {REC['beta']:.6f}, "
            f"`20261012T184537Z_IMP-6-step1-W1b.log:190`); |Im E|/|Re E| = "
            f"{d['imag_ratio']:.3e} ({REC['imag_ratio']:.3e})"
            f"\n[anchor] (ii) Phi_rim,h = {phi_h:.6f}, Phi_rim,cf,h = {phi_cf_h:.6f}, "
            f"|Phi_rim,h/Phi_rim,cf,h - 1| = {phi_miss:.4%} vs MAT4_BAND = "
            f"{m6.MAT4_BAND:.0%} (gate record: {REC['phi_h']} vs {REC['phi_cf_h']}, 1.4181 %)"
            f"\n[anchor] pointwise rim (1/L_x)^2 = {inv_lx**2:.2f} vs the rim-ball mean "
            f"{phi_h:.2f}: averaging over r_r = 2 c^2/a = {2 * c * c / a * 1e3:.2f} mm keeps "
            f"{phi_h * lx**2:.3f} of the peak (predicted 0.358)"
            f"\n[structure] (iii) cavity dofs (reduced) = {d['n_cavity']}, max |E| on them = "
            f"{d['cavity_max']:.3e} (gate record: {REC['n_cavity']}, 0)",
            flush=True,
        )

    assert d["n_cavity"] > 0, "no dofs located on cavity facet tag 2"
    assert d["cavity_max"] < 1e-12, (
        f"max |E| on the cavity dofs is {d['cavity_max']:.3e}, not zero"
    )
    assert beta_miss <= m6.TH15_BETA_BAND, (
        f"beta = {beta:.6f} misses 1 by {beta_miss:.4%} > {m6.TH15_BETA_BAND:.4%} "
        "(imported `TH-15` band) -- an example/test divergence, not a tolerance to move"
    )
    assert phi_miss <= m6.MAT4_BAND, (
        f"Phi_rim,h = {phi_h:.6f} vs Phi_rim,cf,h = {phi_cf_h:.6f}: miss {phi_miss:.4%} "
        f"> {m6.MAT4_BAND:.0%} (imported `MAT-4` band) -- an example/test divergence"
    )

    # ---- printed only: edge against tip, azimuths, rim line -------------------
    big = s >= 2.0 * rho_r
    max_rel_far = float(np.max(np.abs(rel[big]), initial=0.0))
    if comm.rank == 0:
        print(
            f"\n[edge vs tip] same 6 mm long semi-axis, aspect {ASPECT:.0f}: rho_r (disc) = "
            f"{rho_r * 1e3:.3f} mm, rho_t (prolate) = {rho_t * 1e3:.3f} mm; averaging balls "
            f"2c^2/a = {2 * c * c / a * 1e3:.3f} mm, 2a_pro^2/c_pro = "
            f"{2 * a_pro * a_pro / c_pro * 1e3:.3f} mm (predicted 0.24 / 0.24 / 0.48 / 0.48)"
            f"\n[edge vs tip] 1/L_z^pro = {inv_lz_pro:.6f} (record {REC['inv_lz_pro']}) over "
            f"1/L_x^obl = {inv_lx:.6f}: field ratio {inv_lz_pro / inv_lx:.4f} (predicted 2.2350); "
            f"pointwise squares {inv_lz_pro**2:.1f} / {inv_lx**2:.2f}: ratio "
            f"{(inv_lz_pro / inv_lx) ** 2:.3f} (predicted 4.995)"
            f"\n[edge vs tip] ball means: Phi_tip,h record {REC['phi_tip_h']} over this run's "
            f"Phi_rim,h {phi_h:.6f}: ratio {REC['phi_tip_h'] / phi_h:.4f} (predicted 2.1631); "
            f"retained fraction rim {phi_h * lx**2:.3f} (predicted 0.358), tip "
            f"{REC['phi_tip_h'] / inv_lz_pro**2:.3f} (predicted 0.155)",
            flush=True,
        )
        print("[azimuth] |E|/E0 at r = a + rho_r/2, z = 0 (FEM beside closed form):", flush=True)
        for ph, f_, g_ in zip(np.rad2deg(phis), az_fem / E0, az_cf / E0):
            print(
                f"[azimuth] phi = {ph:4.0f} deg: FEM {f_:.4f}, closed form {g_:.4f}, "
                f"rel diff {f_ / g_ - 1.0:+.2%}",
                flush=True,
            )
        mono = bool(np.all(np.diff(az_cf) < 0.0))
        print(
            f"[azimuth] closed form falls monotonically from 0 to 90 deg: {mono} "
            f"(predicted True): two hot arcs at (+-a, 0, 0), none at (0, +-a, 0)",
            flush=True,
        )
        print("[rim line] y = z = 0, x = a + s; |E|/E0 FEM beside closed form (oblate) and the prolate tip line:", flush=True)
        for target in (0.1, 0.5, 1.0, 2.0, 10.0):
            i = int(np.argmin(np.abs(s / rho_r - target)))
            print(
                f"[rim line] s/rho_r = {s[i] / rho_r:7.3f}: FEM {line_fem[i] / E0:8.4f}, "
                f"closed form {line_cf[i] / E0:8.4f}, rel diff {rel[i]:+.2%}, "
                f"prolate tip line {line_pro[i] / E0:8.4f}",
                flush=True,
            )
        print(
            f"[rim line] max |rel diff| for s >= 2 rho_r = {max_rel_far:.2%} (predicted <= 2 %); "
            f"rel diff at s = rho_r/10 = {rel[0]:+.2%} (predicted in [-30 %, +10 %])",
            flush=True,
        )

    # ---- negative control: asserted, same mesh, face-on (axial) drive ----------
    ctl = m6._solve_oblate(ASPECT, "axial", return_fields=True)
    sar_ax, _ = _dg0_sar(msh, ctl.pop("fields"), "sar_ratio_axial")
    delta_beta = abs(ctl["beta"] - 1.0)
    phi_miss_ctl = abs(ctl["phi_h"] / phi_cf_h - 1.0)
    rss_ctl = m6._rss_gib(comm)
    if comm.rank == 0:
        print(
            f"\n[control] axial drive (E0 zhat), same mesh ({ctl['ncells']} cells), "
            f"{ctl['seconds']:.1f} s, RSS(sum) {rss_ctl:.2f} GiB: beta (in-plane basis) = "
            f"{ctl['beta']:.6f} (gate record {REC['ctl_beta']}), |beta - 1| = {delta_beta:.4f} vs "
            f"band {m6.TH15_BETA_BAND:.4%}; beta_axial (axial basis, record, not asserted) = "
            f"{ctl['beta_own']:.6f} ({REC['ctl_beta_axial']})"
            f"\n[control] Phi_rim,h = {ctl['phi_h']:.6f} (record {REC['ctl_phi_h']}) vs the "
            f"in-plane closed form {phi_cf_h:.6f}: miss {phi_miss_ctl:.4%} vs band "
            f"{m6.MAT4_BAND:.0%} (gate record {REC['ctl_miss']} %); cavity dofs "
            f"{ctl['n_cavity']}, max |E| {ctl['cavity_max']:.3e}",
            flush=True,
        )
    assert delta_beta > m6.TH15_BETA_BAND, (
        f"axial beta {ctl['beta']:.6f} is within the band of 1: not a control"
    )
    assert phi_miss_ctl > m6.MAT4_BAND, (
        f"axial Phi miss {phi_miss_ctl:.4%} not over {m6.MAT4_BAND:.0%}: "
        "the Phi gate does not resolve the drive direction"
    )

    # ---- CSV (rank 0), then ParaView (collective) ------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = OUTPUT_DIR / f"{BASENAME}_rim_line.csv"
    if comm.rank == 0:
        header = (
            "s_m,s_over_rho_r,x_m,absE_fem_over_E0,absE_closed_form_over_E0,rel_diff,"
            "Ex_re_fem_over_E0,Ex_closed_form_over_E0,absE_prolate_tip_closed_form_over_E0"
        )
        table = np.column_stack(
            [
                s,
                s / rho_r,
                a + s,
                line_fem / E0,
                line_cf / E0,
                rel,
                line_ex_re / E0,
                line_cf_vec[:, 0] / E0,
                line_pro / E0,
            ]
        )
        np.savetxt(csv_path, table, delimiter=",", header=header, comments="", fmt="%.9e")

    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {
            "sar_ratio_inplane": sar_in,
            "sar_ratio_closed_form": sar_cf_out,
            "sar_ratio_minus_closed_form": sar_diff_out,
            "sar_ratio_axial": sar_ax,
            "rim_mask": rim_mask_out,
            "E_real_inplane": e_real_out,
        },
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[paraview] wrote {xdmf_path}"
            f"\n[paraview] wrote {csv_path}"
            f"\n[paraview] colour by `sar_ratio_inplane` (SAR/SAR0, DG0) on x-z and x-y cuts: the "
            f"hot spots are the two rim arcs at (+-a, 0, 0); `sar_ratio_axial` has none."
            f"\n\nAll assertions hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
