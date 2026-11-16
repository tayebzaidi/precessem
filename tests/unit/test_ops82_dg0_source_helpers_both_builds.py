"""`OPS-82` step 1 (§9 item 849): the `IMP-8` DG0 source helpers build in both
builds.

Brief: docs/planning/chunks/OPS-82.md, "§9 item 849".  ``_rod_source``
(`test_imp8_rod_line_source_fem.py`) and ``_tips_source``
(`test_imp8_rod_tips_fem.py`) now assemble the cell-volume vector as
``one * ufl.conj(w) * ufl.dx``.  On a 12 × 12 × 8 ``create_box`` mesh of
[−0.05, 0.05]² × [−0.03, 0.03] m (no gmsh), in whichever build runs:
(V) the owned DG0 volume sum equals the box volume 6e-4 m³ within 1e-12
relative; (S0) ∫q = P for the rod and P/2 for each ball within ``B_BAND``, the
selected-cell counts > 0.  Controls: complex build — the pre-fix form
``one * w * ufl.dx`` is refused (UFL ``ArityMismatch``, by UFL's own arity
check called directly; see the control's docstring); real build — the
pre-fix form's assembled vector equals the fixed form's bitwise.
"""

from __future__ import annotations

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem, la
from dolfinx import mesh as dmesh
from mpi4py import MPI
from ufl.algorithms.check_arities import ArityMismatch, check_form_arity

from tests.validation.test_imp8_rod_line_source_fem import _rod_source
from tests.validation.test_imp8_rod_tips_fem import _tips_source
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, _allsum, _say

SCALAR_IS_COMPLEX = np.issubdtype(np.dtype(default_scalar_type), np.complexfloating)
BOX_LO = (-0.05, -0.05, -0.03)
BOX_HI = (0.05, 0.05, 0.03)
N_CELLS = (12, 12, 8)
BOX_VOLUME = 0.1 * 0.1 * 0.06  # 6e-4 m³
V_REL_BAND = 1e-12  # the brief's (V) tolerance, relative
POWER = 1.0
LENGTH = 0.03
RADIUS = 0.012


@pytest.fixture(scope="module")
def msh():
    m = dmesh.create_box(MPI.COMM_WORLD, [list(BOX_LO), list(BOX_HI)], list(N_CELLS),
                         cell_type=dmesh.CellType.tetrahedron)
    _say(f"[ops82] build {'complex' if SCALAR_IS_COMPLEX else 'real'}; cells "
         f"{m.topology.index_map(3).size_global}")
    return m


def _dg0(m):
    W = fem.functionspace(m, ("DG", 0))
    w = ufl.TestFunction(W)
    one = fem.Constant(m, np.asarray(1.0, dtype=fem.Function(W).x.array.dtype))
    return W, w, one


def _assembled(form):
    b = fem.assemble_vector(form)
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    return b.array.copy()


def test_v_owned_volume_sum(msh):
    """(V): the helpers' ``bvol`` route, rebuilt from the same fixed form."""
    W, w, one = _dg0(msh)
    arr = _assembled(fem.form(one * ufl.conj(w) * ufl.dx))
    n_loc = msh.topology.index_map(3).size_local
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_loc]
    vol = msh.comm.allreduce(float(np.sum(np.real(arr[dofs]))), op=MPI.SUM)
    rel = abs(vol - BOX_VOLUME) / BOX_VOLUME
    _say(f"[ops82] (V) owned volume sum {vol:.15e} m^3; box {BOX_VOLUME:.15e}; "
         f"rel {rel:.3e} (band {V_REL_BAND:.0e})")
    assert rel <= V_REL_BAND, rel


def test_s0_rod_source(msh):
    q, n_sel, v_sel = _rod_source(msh, power=POWER, length=LENGTH, radius=RADIUS)
    total = np.real(_allsum(msh.comm, q * ufl.dx))
    _say(f"[ops82] (S0) rod: n_sel {n_sel}; v_sel {v_sel:.6e} m^3; int q {total:.15e} W; "
         f"|d| {abs(total - POWER):.3e} (B_BAND {B_BAND:.0e})")
    assert n_sel > 0, n_sel
    assert abs(total - POWER) <= B_BAND, total


def test_s0_tips_source(msh):
    q, balls = _tips_source(msh, power=POWER, radius=RADIUS, length=LENGTH)
    total = np.real(_allsum(msh.comm, q * ufl.dx))
    _say(f"[ops82] (S0) tips: int q {total:.15e} W; |d| {abs(total - POWER):.3e} "
         f"(B_BAND {B_BAND:.0e})")
    assert abs(total - POWER) <= B_BAND, total
    for i, (n_sel, v_sel, qb) in enumerate(balls):
        tb = np.real(_allsum(msh.comm, qb * ufl.dx))
        _say(f"[ops82] (S0) ball {i}: n_sel {n_sel}; v_sel {v_sel:.6e} m^3; int q "
             f"{tb:.15e} W; |d| {abs(tb - 0.5 * POWER):.3e}")
        assert n_sel > 0, (i, n_sel)
        assert abs(tb - 0.5 * POWER) <= B_BAND, (i, tb)


@pytest.mark.skipif(not SCALAR_IS_COMPLEX, reason="complex-build control")
def test_control_complex_prefix_form_refused(msh):
    """The pre-fix form is refused, on every rank, by the UFL arity check that
    ``fem.form`` reaches (`ufl/algorithms/check_arities.py:205`).

    Called directly, not through ``fem.form``: ``ArityMismatch`` derives from
    ``BaseException``, so ``pytest.raises(Exception)`` misses it, and dolfinx's
    ``mpi_jit`` compiles on rank 0 alone, so the raise there left rank 1 hung
    to the 120 s ceiling (`20261106T030126Z_OPS-82-W2-complex.log:119`, Status
    124).  This check is rank-local and needs no JIT."""
    _, w, one = _dg0(msh)
    old = one * w * ufl.dx
    with pytest.raises(ArityMismatch) as exc:
        check_form_arity(old, old.arguments(), complex_mode=True)
    _say(f"[ops82] control (complex): pre-fix form refused, "
         f"{type(exc.value).__name__}: {exc.value}")
    new = one * ufl.conj(w) * ufl.dx
    check_form_arity(new, new.arguments(), complex_mode=True)  # the fixed form passes
    _say("[ops82] control (complex): fixed form passes the same arity check")


@pytest.mark.skipif(SCALAR_IS_COMPLEX, reason="real-build control")
def test_control_real_vector_bitwise_unchanged(msh):
    _, w, one = _dg0(msh)
    old = _assembled(fem.form(one * w * ufl.dx))
    new = _assembled(fem.form(one * ufl.conj(w) * ufl.dx))
    same = bool(np.array_equal(old, new))
    allsame = msh.comm.allreduce(int(same), op=MPI.MIN)
    _say(f"[ops82] control (real): pre-fix vs fixed vector bitwise equal on every rank: "
         f"{bool(allsame)} ({old.size} local entries on rank 0)")
    assert allsame == 1
