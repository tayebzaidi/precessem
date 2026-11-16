"""`WF-13` step 14 (§9 item 731, §10 U113): F-clinical unknowns at any degree.

``_estimate_f_clinical`` returns the first-kind N1curl dimension
``U_k = k E + k(k - 1) F + k(k - 1)(k - 2) T / 2`` from the record's E, F and T, with no
mesh.  Before this step it indexed the record with ``0 if degree == 1 else 1``, so a
degree-3 estimate reported degree-2 unknowns beside NaN prices.

Asserted: k = 1 and k = 2 reproduce the two records exactly, and k = 3 is 22 523 541 from
E, F and T; the old indexing misses U_3 by more than 0.5; on a ``create_box`` mesh, the
identity independent of F-clinical, DolfinX's global ``N1curl`` size equals the formula
at k = 1 … 4.  (NoB) no builder is called.

Not claimed: a degree-3 price or any degree-3 solve.  The prices stay NaN.
"""

from __future__ import annotations

import math

import pytest
from mpi4py import MPI

from tests.package.test_wf13_estimate_f_clinical import _study, no_build  # noqa: F401

E, U2, T = 1_525_193, 7_885_674, 1_147_366
F = (U2 - 2 * E) // 2
U3 = 3 * E + 6 * F + 3 * T


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def test_s14_records_exact_and_degree_3(no_build):  # noqa: F811
    from fem_em_solver import study as study_mod

    assert (E, U2) == study_mod.F_CLINICAL_R1_UNKNOWNS_RECORD
    assert T == study_mod.F_CLINICAL_R1_CELL_RECORD
    assert U2 == 2 * E + 2 * F
    u = {k: _study(degree=k).estimate(MPI.COMM_WORLD).unknowns for k in (1, 2, 3)}
    miss = abs(U2 / u[3] - 1.0)
    _say(f"\n[WF-13 s14] E {E}, F {F}, T {T}; U_k {u}; the old indexing misses U_3 by {miss:.4f}")
    assert u[1] == E and u[2] == U2  # the two records, exactly
    assert u[3] == U3 == 22_523_541
    assert miss > 0.5
    est = _study(degree=3).estimate(MPI.COMM_WORLD)
    assert est.rows_used == () and math.isnan(est.predicted_peak_gib)
    assert math.isnan(est.predicted_seconds)
    assert no_build == []


def test_s14_degree_below_1_refused():
    from fem_em_solver.study import _f_clinical_unknowns

    with pytest.raises(ValueError, match="degree"):
        _f_clinical_unknowns(0)


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_s14_box_identity(k):
    """The formula against DolfinX's own dofmap, on a mesh with nothing clinical in it."""
    from dolfinx import fem, mesh

    from fem_em_solver.study import _n1curl_unknowns

    comm = MPI.COMM_WORLD
    msh = mesh.create_box(comm, [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]], [3, 2, 4],
                          mesh.CellType.tetrahedron)
    for dim in (1, 2):
        msh.topology.create_entities(dim)
    edges, faces, cells = (msh.topology.index_map(d).size_global for d in (1, 2, 3))
    V = fem.functionspace(msh, ("N1curl", k))
    size = V.dofmap.index_map.size_global * V.dofmap.index_map_bs
    formula = _n1curl_unknowns(k, edges, faces, cells)
    _say(f"\n[WF-13 s14] box k={k}: E {edges}, F {faces}, T {cells}; dofmap {size}, formula "
         f"{formula}")
    assert size == formula
