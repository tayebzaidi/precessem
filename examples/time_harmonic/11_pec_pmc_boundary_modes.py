"""Example (`EX-64`): PEC wall against PMC wall on one box.

Every example that names a boundary passes ``"pec_zero_tangential_a"``; `th:8`
and `mri:1` sit on the ``NATURAL`` default, unnamed and unanchored. This one
solves the same box under each option and shows what the wall does: under PEC
the field meets every face normally (``n x E = 0``), under ``NATURAL`` (PMC)
tangentially (``n x curl E = 0``, E itself is free).

**It asserts, and it does not re-implement** (the `ANS-1` rule). The box, the
manufactured fields, the sizes, the bands and the solve are all imported from
``tests/solver/test_boundary_condition_selection.py`` (`TH-3`). The gate's
``_solve_mms`` gains one additive keyword here, ``return_fields``, and its
two literal bands become ``TH3_REL_L2_BAND`` / ``TH3_RATE_FLOOR`` (values
unchanged).

* *Per mode:* owned constrained dofs == ``_th3_exterior_edges(n)`` under PEC,
  == 0 under ``NATURAL``; relative L2 at 24^3 < ``TH3_REL_L2_BAND``, 24^3 below
  12^3, rate > ``TH3_RATE_FLOOR``.
* *Negative control:* E_nat's source under PEC and E_pec's under ``NATURAL``
  at 12^3 each miss by more than ``TH3_REL_L2_BAND``.

Printed, never gated: the gate record beside each reading, and the wall-trace
fraction ``oint |n x E_h|^2 dS / oint |E_h|^2 dS`` (predicted 0 under PEC,
approximately 1 under ``NATURAL``).

Run it::

    ./run_examples.sh -e th:11

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_11_pec_pmc_boundary_modes_pec_combined.xdmf`` and
``..._natural_combined.xdmf`` (CG1 ``E_real``, ``E_error_magnitude``).
"""

from __future__ import annotations

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

from fem_em_solver.core import TimeHarmonicBoundaryCondition  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.solver.test_boundary_condition_selection import (  # noqa: E402
    TH3_RATE_FLOOR,
    TH3_REL_L2_BAND,
    TH3_SIZES,
    _solve_mms,
    _th3_exterior_edges,
    _th3_field_ufl,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_11_pec_pmc_boundary_modes"

PEC = TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A.value
NAT = TimeHarmonicBoundaryCondition.NATURAL.value

#: (key, boundary mode, field key, expected owned constrained dofs, gate record
#: at 12^3 and 24^3 from 20260930T022051Z_TH-3-W1.log:280-289, printed only).
MODES = (
    ("pec", PEC, "pec", _th3_exterior_edges, (1.165799e-01, 5.852485e-02, 0.9942)),
    ("natural", NAT, "nat", lambda n: 0, (7.527560e-02, 3.774071e-02, 0.9961)),
)
#: Control reads at 12^3 (log :293, :297), printed only.
CONTROL_RECORD = {"nat-under-pec": 8.909161e-01, "pec-under-natural": 6.262517e-01}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _reduced(comm, form):
    return comm.allreduce(float(np.real(fem.assemble_scalar(fem.form(form)))), op=MPI.SUM)


def _wall_trace_fraction(e, msh, comm):
    """``oint |n x E|^2 dS / oint |E|^2 dS`` over the box surface, allreduced."""
    n = ufl.FacetNormal(msh)
    tang = _reduced(comm, ufl.inner(ufl.cross(n, e), ufl.cross(n, e)) * ufl.ds)
    full = _reduced(comm, ufl.inner(e, e) * ufl.ds)
    return tang / full


def _export(key, msh, fields, field_key, comm):
    """CG1 interpolants of Re E and |E - E_exact| (an interpolant: nothing asserted)."""
    e = fields.e_complex
    x = ufl.SpatialCoordinate(msh)
    exact = _th3_field_ufl(field_key, x)
    diff = e - exact
    v_vec = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    v_sca = fem.functionspace(msh, ("Lagrange", 1))
    e_real = fem.Function(v_vec, name="E_real")
    e_real.interpolate(fem.Expression(ufl.real(e), v_vec.element.interpolation_points))
    err = fem.Function(v_sca, name="E_error_magnitude")
    err.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(diff, diff))), v_sca.element.interpolation_points)
    )
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_{key}_combined", msh, None,
        {"E_real": e_real, "E_error_magnitude": err}, comm=comm,
    )
    return path


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `th:` group)."
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-64 -- PEC wall against PMC wall: each boundary option, its own field")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] sizes {TH3_SIZES}, band rel L2 < {TH3_REL_L2_BAND}, rate > {TH3_RATE_FLOOR} "
               "(imported from the TH-3 gate module)")

    outputs = []
    figure_done = False
    for key, mode, field_key, expected, record in MODES:
        rel = {}
        for n in TH3_SIZES:
            t = time.perf_counter()
            rel[n], owned, ncells, msh, fields = _solve_mms(mode, field_key, n, return_fields=True)
            _say(comm, f"[{key:7s}] mode={mode} field=E_{field_key} n={n}^3 ({ncells} cells): "
                       f"rel L2 = {rel[n]:.6e}, owned constrained dofs = {owned} "
                       f"(ASSERTED == {expected(n)}) ({time.perf_counter() - t:.1f}s)")
            assert owned == expected(n), (key, n, owned, expected(n))

            if not figure_done:
                # ---- setup figure (EX-57): render on the first mesh built ----
                write_setup_figure(
                    msh, None, FIGURE_DIR / f"{BASENAME}_setup.png",
                    region_names={0: "box (homogeneous gel)"},
                    slice_normal=(0.0, 0.0, 1.0),
                    title="th:11 -- TH-3 box [0, 0.2 m]^3, sigma 0.7 S/m, PEC vs NATURAL (PMC) wall, slice z = L/2",
                    comm=comm,
                )
                figure_done = True

            if n == TH3_SIZES[1]:
                frac = _wall_trace_fraction(fields.e_complex, msh, comm)
                _say(comm, f"[{key:7s}] printed, no band (predicted 0 under PEC, approx 1 under "
                           f"NATURAL): wall-trace fraction oint|n x E_h|^2/oint|E_h|^2 = {frac:.4e}")
                outputs.append(_export(key, msh, fields, field_key, comm))

        coarse, fine = rel[TH3_SIZES[0]], rel[TH3_SIZES[1]]
        rate = float(np.log(coarse / fine) / np.log(TH3_SIZES[1] / TH3_SIZES[0]))
        rc, rf, rr = record
        _say(comm, f"[{key:7s}] rate 12^3 -> 24^3 = {rate:.4f}")
        _say(comm, f"[{key:7s}] gate record (20260930T022051Z_TH-3-W1.log:280-289): {rc:.6e} -> {rf:.6e} "
                   f"({rr:.4f}); this run minus record: {coarse - rc:+.2e}, {fine - rf:+.2e}, "
                   f"{rate - rr:+.4f} (reported, not asserted)")
        assert fine < TH3_REL_L2_BAND, (key, fine)
        assert fine < coarse, (key, coarse, fine)
        assert rate > TH3_RATE_FLOOR, (key, rate)

    # ---- negative control: swapped source / boundary at 12^3 -----------------
    n0 = TH3_SIZES[0]
    for label, mode, field_key in (("nat-under-pec", PEC, "nat"), ("pec-under-natural", NAT, "pec")):
        rel, owned, ncells = _solve_mms(mode, field_key, n0)
        _say(comm, f"[control] {label}: mode={mode} field=E_{field_key} n={n0}^3: rel L2 = {rel:.6e} "
                   f"(ASSERTED > {TH3_REL_L2_BAND}); gate record {CONTROL_RECORD[label]:.6e}")
        assert rel > TH3_REL_L2_BAND, (label, rel)

    for path in outputs:
        _say(comm, f"[output] {path} -- CG1 E_real, E_error_magnitude (interpolants, nothing asserted)")
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
