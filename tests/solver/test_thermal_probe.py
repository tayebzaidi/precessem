"""`MAT-19` step 1 (§9 item 243) — unit checks of the ball-averaged probe.

Box [−1, 1]³, 8³ × 6 tets, P1.  (U1) constant field reads itself; (U2) affine
field reads θ(c) (printed); (U3) a ball containing the box equals ∫θ/∫1 and
volume 8 (rank-safety identity, visible at ``-n 2``); (U3c) planted ghost
control misses (U3) at ``comm.size > 1``; (U4) the reading is bitwise
identical on every rank.
"""

from __future__ import annotations

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.thermal import ProbeBall, ball_average
from fem_em_solver.thermal.probe import _ball_average_on_cells

from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, _allsum


def _box():
    comm = MPI.COMM_WORLD
    msh = dmesh.create_box(comm, [np.array([-1.0, -1.0, -1.0]), np.array([1.0, 1.0, 1.0])],
                           [8, 8, 8], dmesh.CellType.tetrahedron)
    return msh, fem.functionspace(msh, ("Lagrange", 1))


def _fn(V, f):
    u = fem.Function(V)
    u.interpolate(lambda x: f(x).astype(default_scalar_type))
    u.x.scatter_forward()
    return u


def _p(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def test_probe_unit_identities():
    comm = MPI.COMM_WORLD
    msh, V = _box()

    # (U1) constant.
    pb = ProbeBall((0.1, 0.2, 0.05), 0.5)
    u1 = _fn(V, lambda x: np.full(x.shape[1], 3.7))
    v1, vol1 = ball_average(u1, pb)
    nominal = 4.0 * np.pi * 0.5**3 / 3.0
    _p(f"MAT19 (U1) constant 3.7: reading {v1:.15f} rel {abs(v1 / 3.7 - 1):.3e} "
       f"(band {B_BAND:.0e}); V_h/V - 1 = {vol1 / nominal - 1:+.3e} (predicted |.| <= 5e-3)")

    # (U2) affine, printed only.
    u2 = _fn(V, lambda x: 2.0 + x[0] - 0.5 * x[1] + 0.25 * x[2])
    v2, _ = ball_average(u2, pb)
    exact2 = 2.0 + 0.1 - 0.5 * 0.2 + 0.25 * 0.05
    _p(f"MAT19 (U2) affine: reading {v2:.12f} vs theta(c) {exact2:.12f}  rel "
       f"{abs(v2 / exact2 - 1):.3e} (printed; predicted <= 1e-3)")

    # (U3) ball containing the box: the global integral identity.
    u3 = _fn(V, lambda x: 2.0 + x[0] ** 2 + x[1])
    big = ProbeBall((0.0, 0.0, 0.0), 10.0)
    v3, vol3 = ball_average(u3, big)
    one = fem.Constant(msh, default_scalar_type(1.0))
    ref3 = float(np.real(_allsum(comm, u3 * ufl.dx(domain=msh)))) / float(
        np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    e3 = abs(v3 / ref3 - 1.0)
    ev = abs(vol3 / 8.0 - 1.0)
    _p(f"MAT19 (U3) whole-box ball: probe {v3:.15f} vs int/vol {ref3:.15f} rel {e3:.3e}; "
       f"volume {vol3:.15f} rel {ev:.3e} (band {B_BAND:.0e})")

    # (U3c) planted ghost control (collective on every rank).
    e3c = None
    if comm.size > 1:
        imap = msh.topology.index_map(3)
        allc = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
        v3c, vol3c = _ball_average_on_cells(u3, big, allc)
        e3c = max(abs(v3c / ref3 - 1.0), abs(vol3c / 8.0 - 1.0))
        _p(f"MAT19 (U3c) owned+ghost cells: probe rel {abs(v3c / ref3 - 1):.3e}, volume "
           f"{vol3c:.6f} rel {abs(vol3c / 8 - 1):.3e}  (asserted > {B_BAND:.0e}; "
           f"predicted O(1e-2))")

    # (U4) bitwise identical on every rank.
    c0, _ = ball_average(u2, ProbeBall((0.0, 0.0, 0.0), 0.4))
    allv = comm.allgather(c0)
    _p(f"MAT19 (U4) c0 reading on every rank: {allv}")

    assert abs(v1 / 3.7 - 1.0) <= B_BAND, f"(U1) {abs(v1 / 3.7 - 1):.3e}"
    assert e3 <= B_BAND, f"(U3) probe {e3:.3e}"
    assert ev <= B_BAND, f"(U3) volume {ev:.3e}"
    if e3c is not None:
        assert e3c > B_BAND, f"(U3c) ghost-included sum did not miss: {e3c:.3e}"
    assert all(v == allv[0] for v in allv), "(U4) not bitwise identical across ranks"
