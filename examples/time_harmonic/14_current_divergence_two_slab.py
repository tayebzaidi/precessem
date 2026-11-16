"""Example (`EX-71`): does the total current close? ∇·J_tot at CG2 against its vacuous CG1 twin.

Charge conservation outside an impressed source,

    ∇·J_tot = 0,    J_tot = (σ(x) + jωε₀εᵣ(x)) E,

scored by ``fem_em_solver.post.current_divergence`` as a dual norm over the
degree-p Lagrange space vanishing on the wall, normalised by ‖J_tot‖.  The
fixture is `POST-3`'s two-slab box: σ = 0.1 | 1.4 S/m split at x = L/2,
64 MHz, degree-1 N1curl, boundary-driven by the σ_low plane wave.

The point of the example is the pair of test spaces.  ∇(CG1) lies inside the
N1curl₁ test space the solve already enforced its equation against, so the CG1
residual is round-off by Galerkin orthogonality whatever the field is — it
cannot fail.  ∇(CG2) is not inside it, so the CG2 residual is an emergent
property of the discrete solution: first order in h, 9.3e-2 → 6.4e-2 from
n = 8 to n = 12.

Everything gated is **imported** from ``tests/validation/test_current_divergence.py``
(``_solve_piecewise``, ``_residual``, ``COARSE_N``, ``FINE_N`` and the four
bound constants) and ``tests/validation/test_poynting_balance.py``
(``_two_material_mesh`` for the figure, ``_report``, ``OMEGA``, ``MU_R``).
Nothing is restated except the printed records, which the gate holds as log
output, quoted with provenance and never asserted.

Run it through the example runner (the ``th:`` group sources the complex build)::

    ./run_examples.sh -e th:14 -n 2 -t 180
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

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.current_divergence import (  # noqa: E402
    current_divergence_residual,
)
from fem_em_solver.post.power_balance import poynting_power_balance  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_current_divergence import (  # noqa: E402
    CG1_VACUITY_CEILING,
    CG2_CG1_SEPARATION_FLOOR,
    CG2_COARSE_CEILING,
    CG2_RATE_FLOOR,
    COARSE_N,
    FINE_N,
    _residual,
    _solve_piecewise,
)
from tests.validation.test_lossy_plane_wave import BOX_L  # noqa: E402
from tests.validation.test_poynting_balance import (  # noqa: E402
    MU_R,
    OMEGA,
    TAG_HIGH,
    TAG_LOW,
    _report,
    _two_material_mesh,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_14_current_divergence_two_slab"

#: Gate records, 20260930T200032Z_POST-3-closing-W1.log:115,119-124 — printed
#: beside this run, a difference reported and never asserted.
RECORD = {
    "cg2_coarse": 9.316430e-02,
    "cg2_fine": 6.358255e-02,
    "cg1": 6.852517e-15,
    "imbalance_n12": 11.8537,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


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
    _say(comm, "EX-71 -- does the total current close? div J_tot at CG2 vs CG1, two-slab box")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] rate floor {CG2_RATE_FLOOR}, coarse ceiling {CG2_COARSE_CEILING}, "
               f"CG1 ceiling {CG1_VACUITY_CEILING:.0e}, separation floor "
               f"{CG2_CG1_SEPARATION_FLOOR:.0e} (imported from the POST-3 step 3 gate module)")

    # ---- setup figure (EX-57): the figure mesh is its own object -------------
    fig_msh, fig_tags = _two_material_mesh(FINE_N)
    write_setup_figure(
        fig_msh, fig_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={TAG_LOW: "low-σ slab", TAG_HIGH: "high-σ slab"},
        translucent_tags=(TAG_LOW, TAG_HIGH),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.5 * BOX_L, 0.5 * BOX_L, 0.5 * BOX_L),
        title="th:14 -- two-slab box, sigma 0.1 | 1.4 S/m split at x = L/2, "
              "64 MHz; slice z = L/2",
        comm=comm,
    )

    # ---- solves (imported) ----------------------------------------------------
    t = time.perf_counter()
    coarse_solve = _solve_piecewise(COARSE_N)
    fine_solve = _solve_piecewise(FINE_N)
    _say(comm, f"\n[solves] n = {COARSE_N} and n = {FINE_N} ({time.perf_counter() - t:.1f} s)")

    _say(comm, "")
    coarse = _residual(coarse_solve, f"n = {COARSE_N}", degree=2)["relative_residual"]
    fine = _residual(fine_solve, f"n = {FINE_N}", degree=2)["relative_residual"]
    rate = float(np.log(coarse / fine) / np.log(FINE_N / COARSE_N))
    _say(comm, f"\nCG2 rel: coarse {coarse:.6e} (record {RECORD['cg2_coarse']:.6e}), "
               f"fine {fine:.6e} (record {RECORD['cg2_fine']:.6e}); rate in h {rate:.4f} "
               "(reported, not asserted against the record)")

    assert fine < coarse, (coarse, fine)
    _say(comm, f"(1) fine < coarse: ASSERTED ({fine:.4e} < {coarse:.4e})")
    assert rate > CG2_RATE_FLOOR, rate
    _say(comm, f"(2) rate {rate:.4f} > {CG2_RATE_FLOOR}: ASSERTED")
    assert coarse < CG2_COARSE_CEILING, coarse
    _say(comm, f"(3) coarse {coarse:.4e} < {CG2_COARSE_CEILING}: ASSERTED")

    # ---- negative control: CG1 is vacuous, CG2 is not ---------------------------
    _say(comm, "")
    cg1 = _residual(coarse_solve, "vacuity control", degree=1)["relative_residual"]
    cg2 = _residual(coarse_solve, "vacuity control", degree=2)["relative_residual"]
    separation = cg2 / cg1
    _say(comm, f"\n[control] CG1 rel = {cg1:.6e} (gate record {RECORD['cg1']:.6e}, "
               f"round-off, never a record); CG2/CG1 = {separation:.3e} "
               "(gate record 1.36e13)")
    assert cg1 < CG1_VACUITY_CEILING, cg1
    _say(comm, f"[control] CG1 < {CG1_VACUITY_CEILING:.0e}: ASSERTED")
    assert separation > CG2_CG1_SEPARATION_FLOOR, separation
    _say(comm, f"[control] CG2/CG1 > {CG2_CG1_SEPARATION_FLOOR:.0e}: ASSERTED")

    # ---- printed only ------------------------------------------------------------
    fields, sigma_field, eps_field = fine_solve
    msh = sigma_field.function_space.mesh
    assert fields.e_complex.function_space.mesh is msh
    bal = poynting_power_balance(
        fields.e_complex, omega=OMEGA, sigma=sigma_field, mu_r=MU_R, comm=comm
    )
    bal["ncells"] = int(comm.allreduce(msh.topology.index_map(3).size_local, op=MPI.SUM))
    _say(comm, f"\n[printed only] Poynting balance, same n = {FINE_N} field, true sigma(x) "
               f"(record {RECORD['imbalance_n12']} %):")
    if comm.rank == 0:
        _report("true sigma(x)", bal)
    nos = current_divergence_residual(
        coarse_solve[0].e_complex, omega=OMEGA, sigma=coarse_solve[1],
        epsilon_r=coarse_solve[2], degree=2, include_sigma=False, comm=comm,
    )["relative_residual"]
    _say(comm, f"[printed only] include_sigma=False, n = {COARSE_N} CG2: rel = {nos:.6e}, "
               f"ratio to include_sigma=True {nos / coarse:.3f}x "
               "(a non-control by the gate's measurement, about 1.07x)")

    # ---- export: CG1 interpolants (asserted on nothing) + sigma from the solve's mesh ----
    v_vec = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_real = fem.Function(v_vec, name="E_real")
    e_real.interpolate(fem.Expression(ufl.real(fields.e_complex),
                                      v_vec.element.interpolation_points))
    e_imag = fem.Function(v_vec, name="E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(fields.e_complex),
                                      v_vec.element.interpolation_points))
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_n{FINE_N}_combined", msh, None,
        {"E_real": e_real, "E_imag": e_imag, "sigma": sigma_field}, comm=comm,
    )
    _say(comm, f"\n[output] {path}")
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
