"""Steady heat conduction with a volumetric source and a convective boundary.

`MAT-9` step 1 (Phase 8, gel regime).  Solves

.. math::

    -\\nabla\\cdot(k\\nabla T) = Q \\;\\text{in}\\; Ω, \\qquad
    -k\\,\\partial_n T = h\\,(T - T_∞) \\;\\text{on}\\; Γ_R

with Lagrange P1/P2 elements.  The weak form is

.. math::

    a(T, v) = ∫_Ω k∇T·∇v\\,dx + ∮_{Γ_R} h T v\\,ds, \\qquad
    L(v) = ∫_Ω Q v\\,dx + ∮_{Γ_R} h T_∞ v\\,ds .

Because ``v = 1`` lies in every Lagrange space, testing with it gives the
discrete global balance ``∮ h (T_h − T∞) ds = ∫ Q dx`` exactly (to solver
round-off) — the identity the result carries, which checks the Robin sign and
the facet measure rather than the physics.

`MAT-11` adds the Pennes perfusion sink ``w (T − T∞)`` (arterial = ambient,
``w = ω_b ρ_b c_b`` in W/m³/K): ``a`` gains ``∫ w T v dx`` and ``L`` gains
``∫ w T∞ v dx``, only when ``w ≠ 0`` — at ``w = 0`` the forms are `MAT-9`'s.
The ``v = 1`` balance becomes ``∮ h (T_h − T∞) ds + ∫ w (T_h − T∞) dx = ∫ Q dx``.

Real build: the problem is real-valued; nothing here needs complex mode.
Out of scope: EM coupling, T-dependent properties, transient, an arterial
temperature other than the ambient one.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Callable, List, Optional

import dolfinx.fem.petsc as fem_petsc
import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py import PETSc

_PETSC_OPTIONS = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}


@dataclass
class SteadyConductionResult:
    """Solution and its globally reduced balance integrals."""

    temperature: fem.Function
    source_power: float  # ∫ Q dV, allreduced [W]
    convected_power: float  # ∮ h (T − T∞) dS, allreduced [W]
    perfused_power: float = 0.0  # ∫ w (T − T∞) dV, allreduced [W] (`MAT-11`)


def _as_source(mesh, q):
    """A real scalar becomes a ``fem.Constant`` (the `MAT-9` path); a UFL
    expression or ``fem.Function`` is used as-is (`MAT-5` / `MAT-10`)."""
    if isinstance(q, (int, float, np.integer, np.floating)):
        return fem.Constant(mesh, default_scalar_type(q))
    return q


def _solve(
    mesh,
    facet_tags,
    *,
    k: float,
    q,
    h: float,
    t_inf: float,
    robin_tag: int,
    degree: int,
    robin_sign: float = 1.0,
    dirichlet_value: Optional[float] = None,
    petsc_prefix: str = "fem_em_thermal_steady_",
    perfusion_w_per_m3_k: float = 0.0,
) -> SteadyConductionResult:
    """Shared assembly; ``robin_sign`` / ``dirichlet_value`` exist for the
    `MAT-9` negative controls only (a flipped Robin term; a Dirichlet
    ``T = dirichlet_value`` on the Robin surface instead of cooling).  A
    negative ``perfusion_w_per_m3_k`` is accepted here for `MAT-11`'s
    sign-flip control only; the public entry point rejects it."""
    V = fem.functionspace(mesh, ("Lagrange", degree))
    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)

    # A scalar k goes through fem.Constant exactly as before (`MAT-9` forms
    # byte-identical); a UFL expression / DG0 fem.Function is used as-is
    # (`MAT-18` step 2(a), piecewise conductivity).
    k_c = _as_source(mesh, k)
    q_c = _as_source(mesh, q)
    h_c = fem.Constant(mesh, default_scalar_type(h))
    tinf_c = fem.Constant(mesh, default_scalar_type(t_inf))
    s_c = fem.Constant(mesh, default_scalar_type(robin_sign))

    a = k_c * ufl.inner(ufl.grad(T), ufl.grad(v)) * ufl.dx
    L = ufl.inner(q_c, v) * ufl.dx
    bcs = []
    if dirichlet_value is None:
        a = a + ufl.inner(s_c * h_c * T, v) * ds(robin_tag)
        L = L + ufl.inner(s_c * h_c * tinf_c, v) * ds(robin_tag)
    else:
        fdim = mesh.topology.dim - 1
        facets = facet_tags.find(robin_tag)
        dofs = fem.locate_dofs_topological(V, fdim, facets)
        bcs = [fem.dirichletbc(default_scalar_type(dirichlet_value), dofs, V)]
    w_c = None
    if perfusion_w_per_m3_k != 0.0:
        w_c = fem.Constant(mesh, default_scalar_type(perfusion_w_per_m3_k))
        a = a + ufl.inner(w_c * T, v) * ufl.dx
        L = L + ufl.inner(w_c * tinf_c, v) * ufl.dx

    problem = LinearProblem(
        a,
        L,
        bcs=bcs,
        petsc_options=dict(_PETSC_OPTIONS),
        petsc_options_prefix=petsc_prefix,
    )
    T_h = problem.solve()
    if isinstance(T_h, tuple):  # older LinearProblem.solve signatures
        T_h = T_h[0]
    T_h.x.scatter_forward()

    # assemble_scalar is rank-local: reduce every integral.
    comm = mesh.comm
    src = comm.allreduce(fem.assemble_scalar(fem.form(q_c * ufl.dx(domain=mesh))), op=MPI.SUM)
    conv = comm.allreduce(
        fem.assemble_scalar(fem.form(h_c * (T_h - tinf_c) * ds(robin_tag))),
        op=MPI.SUM,
    )
    perf = 0.0
    if w_c is not None:
        perf = comm.allreduce(
            fem.assemble_scalar(fem.form(w_c * (T_h - tinf_c) * ufl.dx(domain=mesh))),
            op=MPI.SUM,
        )
    return SteadyConductionResult(
        T_h, float(np.real(src)), float(np.real(conv)), float(np.real(perf))
    )


def solve_steady_conduction(
    mesh,
    facet_tags,
    *,
    k: float,
    q,
    h: float,
    t_inf: float,
    robin_tag: int,
    degree: int = 1,
    perfusion_w_per_m3_k: float = 0.0,
) -> SteadyConductionResult:
    """Solve ``-∇·(k∇T) + w(T − T∞) = Q`` with ``-k∂ₙT = h(T − T∞)`` on ``robin_tag``.

    Homogeneous ``k``, ``h`` (scalars); ``q`` is a real scalar (wrapped in a
    ``fem.Constant``) or a UFL expression / ``fem.Function`` used as-is in
    ``L`` and in ``source_power``.  Direct solve (MUMPS LU) with
    its own PETSc options prefix.  Returns the temperature and the allreduced
    ``∫Q dV`` and ``∮h(T − T∞) dS``; for a closed Robin boundary the two agree
    to round-off (``v = 1`` test function).

    ``perfusion_w_per_m3_k`` (``w ≥ 0``, default 0: no term is added) is the
    Pennes sink with the arterial temperature equal to ``t_inf``; the result
    then also carries ``perfused_power = ∫w(T − T∞) dV`` and the balance is
    ``convected_power + perfused_power = source_power``.
    """
    if isinstance(k, (int, float, np.integer, np.floating)):
        if k <= 0.0 or h <= 0.0:
            raise ValueError(f"k and h must be positive, got k={k!r}, h={h!r}")
    else:
        if h <= 0.0:
            raise ValueError(f"h must be positive, got h={h!r}")
        if isinstance(k, fem.Function):
            # Rank-local array: allreduce the minimum before checking.
            local = k.x.array
            kmin_loc = float(np.min(np.real(local))) if local.size else np.inf
            kmin = mesh.comm.allreduce(kmin_loc, op=MPI.MIN)
            if not kmin > 0.0:
                raise ValueError(f"k must be positive, got min(k)={kmin!r}")
    if perfusion_w_per_m3_k < 0.0:
        raise ValueError(f"perfusion_w_per_m3_k must be >= 0, got {perfusion_w_per_m3_k!r}")
    return _solve(
        mesh,
        facet_tags,
        k=k,
        q=q,
        h=h,
        t_inf=t_inf,
        robin_tag=robin_tag,
        degree=degree,
        perfusion_w_per_m3_k=perfusion_w_per_m3_k,
    )


_PICARD_PREFIX = itertools.count()


@dataclass
class PicardConductionResult:
    """Picard fixed point of ``-∇·(k∇T) = Q(T)`` and its increment history."""

    temperature: fem.Function
    source_power: float  # ∫ Q(T_{N−1}) dV of the final solve, allreduced [W]
    convected_power: float  # ∮ h (T_N − T∞) dS, allreduced [W]
    increments: List[float] = field(default_factory=list)  # ‖T_{n+1} − T_n‖_L2
    converged: bool = False

    @property
    def iterations(self) -> int:
        return len(self.increments)


def solve_steady_conduction_picard(
    mesh,
    facet_tags,
    *,
    k: float,
    source: Callable,
    h: float,
    t_inf: float,
    robin_tag: int,
    degree: int = 1,
    rtol: float = 1e-12,
    max_iter: int = 60,
    t_initial: Optional[float] = None,
) -> PicardConductionResult:
    """Picard iteration for a temperature-dependent source (`MAT-5` step 1).

    ``source(T)`` maps the current iterate (a ``fem.Function``) to a UFL
    expression Q(T).  Each iteration solves the linear Robin problem with
    ``Q(T_n)`` frozen: the matrix (k-Laplacian + Robin mass) is assembled
    and MUMPS-factored once under its own PETSc options prefix; only ``L`` is
    reassembled.  Starting from ``T_0 = t_initial`` (default ``t_inf``), it
    stops when ``‖T_{n+1} − T_n‖ ≤ rtol·‖T_{n+1} − t_inf‖`` (L2 norms,
    allreduced) or after ``max_iter`` solves.  For ``Q = q₀[1 + α(T − T∞)]``
    the increments contract by ``(mR/λ₁)²`` on the Robin ball.

    Returns the final iterate, the allreduced balance terms of the final
    solve (``∫Q(T_{N−1})`` vs ``∮h(T_N − T∞)`` — equal to round-off by the
    ``v = 1`` test function) and the increment history.
    """
    if k <= 0.0 or h <= 0.0:
        raise ValueError(f"k and h must be positive, got k={k!r}, h={h!r}")
    comm = mesh.comm
    V = fem.functionspace(mesh, ("Lagrange", degree))
    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)

    k_c = fem.Constant(mesh, default_scalar_type(k))
    h_c = fem.Constant(mesh, default_scalar_type(h))
    tinf_c = fem.Constant(mesh, default_scalar_type(t_inf))

    T_old = fem.Function(V)
    T_new = fem.Function(V)
    T_old.x.array[:] = t_inf if t_initial is None else t_initial
    T_old.x.scatter_forward()
    q_expr = _as_source(mesh, source(T_old))

    a_form = fem.form(
        k_c * ufl.inner(ufl.grad(T), ufl.grad(v)) * ufl.dx + ufl.inner(h_c * T, v) * ds(robin_tag)
    )
    L_form = fem.form(ufl.inner(q_expr, v) * ufl.dx + ufl.inner(h_c * tinf_c, v) * ds(robin_tag))
    src_form = fem.form(q_expr * ufl.dx(domain=mesh))
    conv_form = fem.form(h_c * (T_new - tinf_c) * ds(robin_tag))
    inc_form = fem.form((T_new - T_old) ** 2 * ufl.dx)
    nrm_form = fem.form((T_new - tinf_c) ** 2 * ufl.dx)

    def _allsum(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    A = fem_petsc.assemble_matrix(a_form)
    A.assemble()
    b = T_new.x.petsc_vec.duplicate()
    prefix = f"fem_em_thermal_picard_{next(_PICARD_PREFIX)}_"
    ksp = PETSc.KSP().create(comm)
    ksp.setOptionsPrefix(prefix)
    opts = PETSc.Options()
    for key, val in _PETSC_OPTIONS.items():
        opts[prefix + key] = val
    ksp.setFromOptions()
    ksp.setOperators(A)

    result = PicardConductionResult(T_new, 0.0, 0.0)
    try:
        for _ in range(int(max_iter)):
            with b.localForm() as b_local:
                b_local.set(0.0)
            fem_petsc.assemble_vector(b, L_form)
            b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            ksp.solve(b, T_new.x.petsc_vec)
            T_new.x.scatter_forward()
            inc = np.sqrt(_allsum(inc_form))
            nrm = np.sqrt(_allsum(nrm_form))
            result.increments.append(float(inc))
            # Balance of this solve: the source is still Q(T_old).
            result.source_power = _allsum(src_form)
            result.convected_power = _allsum(conv_form)
            T_old.x.array[:] = T_new.x.array
            T_old.x.scatter_forward()
            if inc <= rtol * nrm:
                result.converged = True
                break
    finally:
        ksp.destroy()
        A.destroy()
        b.destroy()
        for key in _PETSC_OPTIONS:
            opts.delValue(prefix + key)
    return result
