"""Example (`EX-82`): the cracked TE101 box — two Leontovich lips against Pozar.

`th:10` (`EX-60`) is the uncracked lossy box. Here the y = 0.4 plane of the
same 1.0 x 0.8 x 0.6 m box is a **crack**: two separate half-cavities, each lip
an exterior facet with its own outward normal, so `TH-14`'s ``ds`` Leontovich
term acts on both sides. Beside it, the same plane as a *single-valued* sheet
(the merged mesh with ``sheet=(septum, SEPTUM_TAG)``) is the negative control:
TE101 has E_t = 0 on y = 0.4, so that sheet cannot see the mode.

**Imported, never restated** (the `ANS-1` rule): every constant, band and
helper comes from ``tests/validation/test_geo33_two_sided_leontovich.py``
(`GEO-33` step 0b), Pozar's ``pozar_te10l_q_conductor`` and
``core.cracked_box``. The solver is the gate's own
``solve_impedance_mesh_mode``, which gains one additive keyword here,
``return_vectors`` (and the defaulted ``MeshImpedanceModes.mode_functions``),
so the example can export the eigenfunctions whose Q it asserts.

Asserted at the imported bands: (Q0) merged walls vs Pozar; (Qs) both members
of the cracked pair vs Pozar of the half box; (R) Q_sep/Q0 vs ``RATIO_LITERAL``;
(F) and (D) on every solve; the exported-mode identity (the Rayleigh quotient
returns λ to 1e-6 relative for the Q0 mode, both pair members and the
control); the negative control. Printed, no band: the lip share of the wall
|n x E|^2, the septum E_t fraction, the lower-half energy fraction of each pair
member (basis dependent, not a label) and 1/Q_sep - 1/Q0.

Run it::

    ./run_examples.sh -e th:15

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_15_cracked_box_two_sided_lips_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.core.cavity import (  # noqa: E402
    pozar_te10l_q_conductor,
    solve_impedance_mesh_mode,
)
from fem_em_solver.core.cracked_box import (  # noqa: E402
    LIP_LOWER,
    LIP_UPPER,
    Y_PLANE,
    cracked_box_mesh,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.constants import C_0, MU_0  # noqa: E402

from tests.validation import test_cavity_leontovich_q, test_cavity_resonances  # noqa: E402
from tests.validation import test_geo33_two_sided_leontovich as g33  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_15_cracked_box_two_sided_lips"

RAYLEIGH_RTOL = 1e-6
LIP_SHARE_PREDICTED = (0.2274, 0.2265, 0.2290)  # centre, low, high
SEPTUM_FRACTION_PREDICTED = 1e-8
TOL_BAND = g33.Q_TOLERANCE


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _red(comm, form):
    return comm.allreduce(complex(fem.assemble_scalar(fem.form(form))), op=MPI.SUM)


def _nxe_sq(mode, n, measure):
    return ufl.inner(ufl.cross(n, mode), ufl.cross(n, mode)) * measure


def _rayleigh(mode, msh, z_s, comm, ctags=None, septum=None):
    """(∫|curl E|² + γ∫_ds|n×E|² [+ γΣ±∫_dS(septum)|n±×E±|²]) / ∫|E|², reduced."""
    omega0 = C_0 * np.sqrt(g33.K0_SQ)
    gamma = fem.Constant(msh, default_scalar_type(1j * omega0 * MU_0 / complex(z_s)))
    n = ufl.FacetNormal(msh)
    num = _red(comm, ufl.inner(ufl.curl(mode), ufl.curl(mode)) * ufl.dx)
    num += _red(comm, gamma * _nxe_sq(mode, n, ufl.ds))
    if septum is not None:
        dS = ufl.Measure("dS", domain=msh, subdomain_data=septum)
        for s in ("+", "-"):
            num += _red(comm, gamma * ufl.inner(
                ufl.cross(n(s), mode(s)), ufl.cross(n(s), mode(s))) * dS(g33.SEPTUM_TAG))
    return num / _red(comm, ufl.inner(mode, mode) * ufl.dx)


def _dg0_sq(msh, expr_sq, name):
    dg0 = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(dg0, name=name)
    f.interpolate(fem.Expression(expr_sq, dg0.element.interpolation_points))
    return dg0, f


def _normalised_magnitude(msh, expr_sq, name, comm):
    """sqrt of a DG0-interpolated ``inner(.,.)``, in numpy, scaled by its allreduced max."""
    dg0, f = _dg0_sq(msh, expr_sq, name)
    mag = np.sqrt(np.maximum(np.real(f.x.array), 0.0))
    n_owned = dg0.dofmap.index_map.size_local
    peak = comm.allreduce(float(np.max(mag[:n_owned])) if n_owned else 0.0, op=MPI.MAX)
    out = fem.Function(dg0, name=name)
    out.x.array[:] = mag / peak
    out.x.scatter_forward()
    return out


def _half_index(msh):
    dg0 = fem.functionspace(msh, ("DG", 0))
    n_cells = msh.topology.index_map(3).size_local + msh.topology.index_map(3).num_ghosts
    mid = dmesh.compute_midpoints(msh, 3, np.arange(n_cells, dtype=np.int32))
    cell_to_dof = np.asarray(dg0.dofmap.list).reshape(-1)
    h = fem.Function(dg0, name="half_index")
    h.x.array[cell_to_dof] = (mid[:, 1] > Y_PLANE).astype(np.float64)
    h.x.scatter_forward()
    return h


def _check(comm, label, res, i=0):
    q = g33._check_mode(label, res, i)  # prints via the gate's own _say (rank 0)
    return q


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this for `th:`)."
        )
    _say(comm, "=" * 78)
    _say(comm, "EX-82 -- the cracked TE101 box: two-sided Leontovich lips against Pozar")
    _say(comm, "=" * 78)

    # Structural identities (printed): the bands are the gate modules' own objects.
    is_q = g33.Q_TOLERANCE is test_cavity_leontovich_q.Q_TOLERANCE
    is_f = g33.FREQUENCY_TOLERANCE_PCT is test_cavity_resonances.FREQUENCY_TOLERANCE_PCT
    _say(comm, f"[imports] g33.Q_TOLERANCE is test_cavity_leontovich_q.Q_TOLERANCE: {is_q}; "
               f"g33.FREQUENCY_TOLERANCE_PCT is test_cavity_resonances.FREQUENCY_TOLERANCE_PCT: "
               f"{is_f}")
    assert is_q and is_f

    merged, _mtags, spread_m, septum = cracked_box_mesh(
        merge=True, comm=comm, septum_tag=g33.SEPTUM_TAG)
    cracked, ctags, spread_c = cracked_box_mesh(merge=False, comm=comm)
    assert merged is not None and cracked is not None, (spread_m, spread_c)

    # ---- setup figure (EX-57), right after the mesh is built ----
    write_setup_figure(
        cracked, None, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={0: "cavity (vacuum), cracked at y = 0.4"},
        slice_normal=(0.0, 0.0, 1.0),
        show_edges_on_slice=True,
        title="th:15 -- TE101 box cracked at y = 0.4: two half-cavities, two lips",
        comm=comm,
    )

    lo = g33._area(cracked, ctags, LIP_LOWER)
    up = g33._area(cracked, ctags, LIP_UPPER)
    _say(comm, f"[lips] areas lower={lo!r} upper={up!r} beside 0.6 (the gate's :90)")

    q0_pozar = pozar_te10l_q_conductor(g33.EDGES, 1, g33.SIGMA_GATE)
    qs_pozar = pozar_te10l_q_conductor(g33.HALF_EDGES, 1, g33.SIGMA_GATE)

    # ---- (Q0) merged walls only ----
    r0 = solve_impedance_mesh_mode(
        merged, sigma_s_per_m=g33.SIGMA_GATE, k0_sq=g33.K0_SQ, degree=g33.DEGREE,
        return_vectors=True)
    q0 = _check(comm, "Q0", r0)
    e0 = abs(q0 / q0_pozar - 1.0)
    _say(comm, f"[Q0] Q0={q0:.6e} vs Pozar {q0_pozar:.6e}: {100 * e0:.4f} % "
               f"(record 8.017516e+02 vs 8.016849e+02, 0.0083 %; ASSERTED < {100 * TOL_BAND:g} %)")

    # ---- (Qs) cracked, both members of the pair ----
    rs = solve_impedance_mesh_mode(
        cracked, sigma_s_per_m=g33.SIGMA_GATE, k0_sq=g33.K0_SQ, degree=g33.DEGREE,
        return_vectors=True)
    assert rs.n_converged >= 2, rs.n_converged
    qs_a = _check(comm, "Qs", rs, 0)
    qs_b = _check(comm, "Qs", rs, 1)
    es_a, es_b = abs(qs_a / qs_pozar - 1.0), abs(qs_b / qs_pozar - 1.0)
    _say(comm, f"[Qs] members {qs_a:.6e} / {qs_b:.6e} vs Pozar {qs_pozar:.6e}: "
               f"{100 * es_a:.4f} % / {100 * es_b:.4f} % (record 6.192100e+02, 0.0611 %)")
    ratio = qs_a / q0
    er = abs(ratio / g33.RATIO_LITERAL - 1.0)
    _say(comm, f"[R] Q_sep/Q0={ratio:.6f} vs {g33.RATIO_LITERAL}: {100 * er:.4f} % "
               f"(record 0.772322, 0.0528 %)")
    dinv = 1.0 / qs_a - 1.0 / q0
    _say(comm, f"[printed] 1/Q_sep - 1/Q0 = {dinv:.6e} beside 3.676919e-04")

    # ---- negative control: merged mesh + single-valued septum sheet ----
    rc = solve_impedance_mesh_mode(
        merged, sigma_s_per_m=g33.SIGMA_GATE, k0_sq=g33.K0_SQ, degree=g33.DEGREE,
        sheet=(septum, g33.SEPTUM_TAG), return_vectors=True)
    qc = _check(comm, "control", rc)
    rc_ratio = qc / q0
    miss = abs(rc_ratio / g33.RATIO_LITERAL - 1.0)
    ceiling = 1.0 / g33.RATIO_LITERAL - 1.0
    _say(comm, f"[control] Q_sheet/Q0={rc_ratio:.8f} (|delta| {abs(rc_ratio - 1):.3e}); miss from "
               f"{g33.RATIO_LITERAL}: {100 * miss:.3f} % (ceiling {100 * ceiling:.3f} %, "
               f"{miss / TOL_BAND:.2f}x band)  [ASSERTED: within band of 1, outside "
               f"(R)'s band, <= ceiling]")

    assert e0 < g33.Q_TOLERANCE, (q0, q0_pozar)
    assert es_a < g33.Q_TOLERANCE and es_b < g33.Q_TOLERANCE, (qs_a, qs_b, qs_pozar)
    assert er < g33.Q_TOLERANCE, ratio
    assert abs(rc_ratio - 1.0) < g33.Q_TOLERANCE, rc_ratio
    assert miss > g33.Q_TOLERANCE, rc_ratio
    assert miss <= ceiling * (1.0 + 1e-6), (miss, ceiling)

    # ---- the exported-mode identity ----
    cases = (
        ("Q0 mode", r0.mode_functions[0], merged, r0.eigenvalues[0], r0.surface_impedance_ohm, None),
        ("pair a", rs.mode_functions[0], cracked, rs.eigenvalues[0], rs.surface_impedance_ohm, None),
        ("pair b", rs.mode_functions[1], cracked, rs.eigenvalues[1], rs.surface_impedance_ohm, None),
        ("control", rc.mode_functions[0], merged, rc.eigenvalues[0], rc.surface_impedance_ohm, septum),
    )
    rq_worst = 0.0
    for label, mode, msh, lam, z_s, sep in cases:
        lam_rq = _rayleigh(mode, msh, z_s, comm, septum=sep)
        rel = abs(lam_rq - lam) / abs(lam)
        rq_worst = max(rq_worst, rel)
        _say(comm, f"[mode] {label}: Rayleigh {lam_rq:.10e} vs lambda {complex(lam):.10e}: "
                   f"rel {rel:.3e} (ASSERTED <= {RAYLEIGH_RTOL:g})")
        assert rel <= RAYLEIGH_RTOL, (label, lam_rq, lam, rel)

    # ---- printed: lip share, septum fraction, half energy fraction ----
    n_c = ufl.FacetNormal(cracked)
    ds_c = ufl.Measure("ds", domain=cracked, subdomain_data=ctags)
    h = _half_index(cracked)
    c_lo, c_hi, c_ctr = (LIP_SHARE_PREDICTED[1], LIP_SHARE_PREDICTED[2], LIP_SHARE_PREDICTED[0])
    for label, mode in (("pair a", rs.mode_functions[0]), ("pair b", rs.mode_functions[1])):
        lips = _red(comm, _nxe_sq(mode, n_c, ds_c(LIP_LOWER)) + _nxe_sq(mode, n_c, ds_c(LIP_UPPER))).real
        total = _red(comm, _nxe_sq(mode, n_c, ds_c)).real
        share = lips / total
        verdict = "inside" if c_lo <= share <= c_hi else "OUTSIDE (reported, not a stop)"
        e_tot = _red(comm, ufl.inner(mode, mode) * ufl.dx).real
        e_up = _red(comm, ufl.inner(mode, mode) * h * ufl.dx).real
        _say(comm, f"[printed] lip share s_h ({label}) = {share:.6f}; predicted {c_ctr} "
                   f"[{c_lo}, {c_hi}] -> {verdict}; beside 1 - R(measured) = {1 - ratio:.6f}, "
                   f"1 - RATIO_LITERAL = {1 - g33.RATIO_LITERAL:.6f}")
        _say(comm, f"[printed] lower-half energy fraction ({label}) = {1 - e_up / e_tot:.6f} "
                   "(basis dependent; not a label)")
    n_m = ufl.FacetNormal(merged)
    dS_m = ufl.Measure("dS", domain=merged, subdomain_data=septum)
    cm = rc.mode_functions[0]
    sep_et = sum(
        _red(comm, ufl.inner(ufl.cross(n_m(s), cm(s)), ufl.cross(n_m(s), cm(s))) * dS_m(g33.SEPTUM_TAG)).real
        for s in ("+", "-"))
    wall_et = _red(comm, _nxe_sq(cm, n_m, ufl.ds)).real
    _say(comm, f"[printed] septum E_t fraction (control mode) = {sep_et / wall_et:.3e} "
               f"(predicted <= {SEPTUM_FRACTION_PREDICTED:g})")

    # ---- ParaView ----
    OUTPUT_DIR.mkdir(exist_ok=True)
    ea = _normalised_magnitude(cracked, ufl.inner(rs.mode_functions[0], rs.mode_functions[0]),
                               "E_magnitude_pair_a", comm)
    eb = _normalised_magnitude(cracked, ufl.inner(rs.mode_functions[1], rs.mode_functions[1]),
                               "E_magnitude_pair_b", comm)
    ca = _normalised_magnitude(
        cracked, ufl.inner(ufl.curl(rs.mode_functions[0]), ufl.curl(rs.mode_functions[0])),
        "curlE_magnitude_pair_a", comm)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", cracked, None,
        {"E_magnitude_pair_a": ea, "E_magnitude_pair_b": eb,
         "curlE_magnitude_pair_a": ca, "half_index": h},
        comm=comm, facet_tags=ctags)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- E_magnitude_pair_a/_b, curlE_magnitude_pair_a, half_index "
               "(DG0; facet block: lips 11/12, faces 1-6)")
    _say(comm, f"[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
