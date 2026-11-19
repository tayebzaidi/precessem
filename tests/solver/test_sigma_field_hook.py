"""`MAT-5` step 3 (§9 item 237) — ``TimeHarmonicSolver.solve(sigma_field_hook=…)``.

The knob is additive and unset-identical.  On ``test_time_harmonic_smoke``'s
fixture (``_smoke_mesh(0.03)``, SIGMA 0.7, 127.74 MHz, axial drive on tag 1):

(a) no hook against an identity hook (each owned+ghost cell rewritten with its
    own value): the owned E coefficients are ``np.array_equal`` (in-process,
    defined at ``-n 1``, `OPS-43`);
(b) a doubling hook: the returned ``fields.sigma_field`` reads exactly
    2 × 0.7 on every owned cell, i.e. the solved field is the hooked one;
(c) control (N3): the doubling hook moves ‖E‖ by more than 1e-6 relative
    (predicted O(1e-1), printed).
"""

from __future__ import annotations

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem, TimeHarmonicSolver
from tests.complex_mode import complex_only
from tests.solver.test_time_harmonic_smoke import (
    EPSILON_R,
    FREQUENCY_HZ,
    SIGMA,
    _smoke_mesh,
)

N3_FLOOR = 1.0e-6


def _solve(hook=None):
    comm = MPI.COMM_WORLD
    mesh, cell_tags, facet_tags = _smoke_mesh(0.03, comm)
    problem = TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA, epsilon_r=EPSILON_R, mu_r=1.0),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
    )
    solver = TimeHarmonicSolver(problem, degree=1)
    return solver.solve(
        current_density=lambda x: ufl.as_vector([0.0, 0.0, 1.0]),
        subdomain_id=1,
        gauge_penalty=1e-3,
        sigma_field_hook=hook,
    )


def _owned(f):
    imap = f.function_space.dofmap.index_map
    n = imap.size_local * f.function_space.dofmap.index_map_bs
    return f.x.array[:n]


def _norm(f):
    comm = f.function_space.mesh.comm
    return float(np.sqrt(comm.allreduce(float(np.sum(np.abs(_owned(f)) ** 2)), op=MPI.SUM)))


def _identity(sf):
    sf.x.array[:] = sf.x.array.copy()
    sf.x.scatter_forward()


def _double(sf):
    sf.x.array[:] = 2.0 * sf.x.array
    sf.x.scatter_forward()


@pytest.fixture(scope="module")
def solved():
    return {"none": _solve(None), "ident": _solve(_identity), "double": _solve(_double)}


@complex_only
def test_a_identity_hook_is_bitwise_unset(solved):
    a = _owned(solved["none"].e_complex)
    b = _owned(solved["ident"].e_complex)
    same = bool(np.array_equal(a, b))
    comm = MPI.COMM_WORLD
    same = comm.allreduce(int(same), op=MPI.MIN) == 1
    if comm.rank == 0:
        print(
            f"\nMAT5s3 W0 (a) identity hook vs no hook: array_equal = {same}  "
            f"(owned dofs on rank 0: {a.size})",
            flush=True,
        )
    assert same, "(a) identity hook changed the owned E coefficients"


@complex_only
def test_b_doubling_hook_is_the_solved_sigma(solved):
    sf = solved["double"].sigma_field
    n = sf.function_space.dofmap.index_map.size_local
    loc = float(np.max(np.abs(sf.x.array[:n] - 2.0 * SIGMA))) if n else 0.0
    err = MPI.COMM_WORLD.allreduce(loc, op=MPI.MAX)
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\nMAT5s3 W0 (b) max|sigma_field - 2*{SIGMA}| over owned cells = {err:.3e} "
            f"(exact 0 required)",
            flush=True,
        )
    assert err == 0.0, f"(b) fields.sigma_field is not the hooked field: {err:.3e}"


@complex_only
def test_c_doubling_hook_moves_e(solved):
    n0 = _norm(solved["none"].e_complex)
    diff = solved["double"].e_complex.x.array - solved["none"].e_complex.x.array
    comm = MPI.COMM_WORLD
    nd = float(
        np.sqrt(
            comm.allreduce(
                float(np.sum(np.abs(diff[: _owned(solved["none"].e_complex).size]) ** 2)),
                op=MPI.SUM,
            )
        )
    )
    rel = nd / n0
    if comm.rank == 0:
        print(
            f"\nMAT5s3 W0 (c)/(N3) ||E_2sigma - E||/||E|| = {rel:.6e}  "
            f"(asserted > {N3_FLOOR:.0e}, predicted O(1e-1))",
            flush=True,
        )
    assert rel > N3_FLOOR, f"(N3) doubling hook moved E by only {rel:.3e}"
