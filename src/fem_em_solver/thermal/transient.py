"""Transient heat conduction with a volumetric source and a convective boundary.

`MAT-9` step 2 (Phase 8, gel regime).  Solves

.. math::

    ρc\\,∂_t T - \\nabla\\cdot(k\\nabla T) = Q \\;\\text{in}\\; Ω, \\qquad
    -k\\,\\partial_n T = h\\,(T - T_∞) \\;\\text{on}\\; Γ_R,

with Lagrange elements in space and the θ-scheme in time (θ = 1 backward
Euler, θ = ½ Crank–Nicolson; Jin, *The FEM in Electromagnetics* 3e,
§12.2.2–12.3 for the first order and stability of BE and the θ family):

.. math::

    \\Big(\\frac{ρc}{Δt}M + θ(K + H)\\Big)T^{n+1}
    = \\Big(\\frac{ρc}{Δt}M - (1-θ)(K + H)\\Big)T^{n} + F,

with ``F = ∫Q v + ∮h T∞ v`` (step source: ``Q`` is on for all ``t > 0``).
The system matrix is assembled and MUMPS-factored once per ``(Δt, θ)`` under
its own PETSc prefix; each step rebuilds only the right-hand side (the
held-factor pattern of ``core/time_harmonic.py``).  ``startup_half_steps``
(even) replaces the first ``startup_half_steps/2`` full steps by backward
Euler half-steps — Rannacher's start, which damps the stiff modes CN leaves
undamped.

Testing with ``v = 1`` (in every Lagrange space; ``K·1 = 0``) gives the
per-step discrete balance

    ρc∫(Tⁿ⁺¹ − Tⁿ)/Δt dx + ∮h(T^(n+θ) − T∞) ds − ∫Q dx = 0

to solver round-off; its worst relative value over all steps is returned.

`MAT-12`: ``q`` may be a UFL expression / ``fem.Function`` (e.g. a DG0 SAR
source), and every test-function product is written ``ufl.inner(·, v)`` so the
forms also compile in the complex build (the problem itself stays real).
`MAT-13`: an optional Pennes sink ``w(T − T∞)`` (``perfusion_w_per_m3_k``),
θ-weighted with the Robin term.  `MAT-17` step 2: an optional per-step source
schedule ``χ(t)`` (``source_schedule``) multiplies the source term; it is
evaluated at the midpoint of every (half-)step, with step times taken from
the integer step index.  With ``source_schedule=None`` (the default) every
form is the unscheduled one.  Out of scope: T-dependent properties,
adaptive Δt.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Sequence

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
import dolfinx.fem.petsc as fem_petsc
from mpi4py import MPI
from petsc4py import PETSc

from .steady import _as_source

_PREFIX_COUNTER = itertools.count()


def _is_scalar(x) -> bool:
    return isinstance(x, (int, float, np.integer, np.floating))


def _check_positive_coefficient(mesh, name, f) -> None:
    """A ``fem.Function`` coefficient (`MAT-18` step 2(b), DG0 k / ρc) must be
    positive: the real part's minimum is allreduced (the array is rank-local).
    Other UFL expressions are not checked."""
    if isinstance(f, fem.Function):
        local = f.x.array
        fmin_loc = float(np.min(np.real(local))) if local.size else np.inf
        fmin = mesh.comm.allreduce(fmin_loc, op=MPI.MIN)
        if not fmin > 0.0:
            raise ValueError(f"{name} must be positive, got min({name})={fmin!r}")


@dataclass
class TransientConductionResult:
    """Final temperature, requested snapshots and the balance record."""

    temperature: fem.Function
    snapshots: Dict[int, fem.Function] = field(default_factory=dict)
    max_balance_residual: float = 0.0
    balance_residuals: list = field(default_factory=list)
    source_power: float = 0.0  # ∫ Q dV, allreduced [W]
    n_factorizations: int = 0
    perfused_power: float = 0.0  # ∫ w (T^(n+θ) − T∞) dV of the last step [W]
    source_multipliers: list = field(default_factory=list)  # χₙ per (half-)step; schedule only


def solve_transient_conduction(
    mesh,
    facet_tags,
    *,
    k,
    rho_c,
    q,
    h: float,
    t_inf: float,
    t_init: float,
    robin_tag: int,
    dt: float,
    n_steps: int,
    theta: float = 1.0,
    startup_half_steps: int = 0,
    degree: int = 1,
    snapshot_steps: Sequence[int] = (),
    perfusion_w_per_m3_k: float = 0.0,
    source_schedule: Optional[Callable[[float], float]] = None,
    t_init_expr: Optional[Callable] = None,
) -> TransientConductionResult:
    """θ-scheme for ``ρc ∂ₜT − ∇·(k∇T) + w(T − T∞) = Q`` with Robin cooling on
    ``robin_tag``.

    ``q`` is a real scalar (wrapped in a ``fem.Constant``) or a UFL expression
    / ``fem.Function`` used as-is (`steady._as_source`).  So are ``k`` and
    ``rho_c`` (`MAT-18` step 2(b): piecewise DG0 coefficients); a scalar takes
    the unchanged ``fem.Constant`` path, a ``fem.Function`` is checked
    positive by an allreduced minimum, and the same ``rho_c`` object weights
    the system matrix, the right-hand side and the balance's storage term.
    ``n_steps`` full steps of ``dt`` from ``T(0) = t_init`` (uniform).  With
    ``startup_half_steps = m`` (even) the first ``m/2`` full steps are taken as
    ``m`` backward-Euler half-steps.  ``snapshot_steps`` are full-step indices
    (``0`` is the initial state); snapshots are copies of ``T`` there.  The
    balance residual of every (half-)step is
    ``|ρc∫(Tⁿ⁺¹−Tⁿ)/Δt + ∮h(T^(n+θ)−T∞) + ∫w(T^(n+θ)−T∞) − ∫Q| / ∫Q``, all
    integrals allreduced.

    ``perfusion_w_per_m3_k`` (``w ≥ 0``, default 0: no term is added, the
    forms are `MAT-9`'s) is the Pennes sink with the arterial temperature equal
    to ``t_inf`` (`MAT-13`), θ-weighted on both sides like the Robin term; the
    result then carries the last step's ``perfused_power``.

    ``source_schedule`` (default ``None``: the source is on for all ``t > 0``
    and the forms are unchanged) is a callable ``χ(t_mid) -> float`` (`MAT-17`
    step 2) evaluated at the midpoint of every (half-)step: half-step ``i``
    spans ``[iΔt/2, (i+1)Δt/2]`` and full step ``i`` after the start spans
    ``[(m/2+i)Δt, (m/2+i+1)Δt]``.  It multiplies ``∫Q v`` for that step.  The
    balance residual is then ``|stor + conv + perf − χₙ∫Q| / ∫Q`` (normalised
    by the full-on ``source_power``, never by ``χₙ∫Q``, which vanishes on an
    off step), and ``source_multipliers`` records ``χₙ`` per (half-)step.

    ``t_init_expr`` (keyword-only, default ``None``: ``T(0) = t_init``,
    unchanged) is a callable ``x -> T(x)`` (DolfinX ``interpolate`` signature,
    ``x`` of shape ``(3, n)``) interpolated into the Lagrange space as the
    initial state; when given it overrides ``t_init`` (`MAT-21` step 19,
    §9 item 669).  When the source integral ``∫Q`` is exactly zero (where the
    old normalisation divides by zero) the per-step residual is instead
    normalised by the initial stored energy per step,
    ``∫ρc|T(0) − T∞| dV / Δtₙ`` (allreduced; the same units as the residual),
    so ``|stor + conv + perf| · Δtₙ / ∫ρc|T(0) − T∞|``.  With ``∫Q ≠ 0`` the
    old normalisation is used unchanged.
    """
    if perfusion_w_per_m3_k < 0.0:
        raise ValueError(
            f"perfusion_w_per_m3_k must be >= 0, got {perfusion_w_per_m3_k!r}"
        )
    return _solve_transient(
        mesh, facet_tags, k=k, rho_c=rho_c, q=q, h=h, t_inf=t_inf, t_init=t_init,
        robin_tag=robin_tag, dt=dt, n_steps=n_steps, theta=theta,
        startup_half_steps=startup_half_steps, degree=degree,
        snapshot_steps=snapshot_steps, perfusion_w_per_m3_k=perfusion_w_per_m3_k,
        source_schedule=source_schedule, t_init_expr=t_init_expr,
    )


def _solve_transient(
    mesh,
    facet_tags,
    *,
    k,
    rho_c,
    q,
    h: float,
    t_inf: float,
    t_init: float,
    robin_tag: int,
    dt: float,
    n_steps: int,
    theta: float = 1.0,
    startup_half_steps: int = 0,
    degree: int = 1,
    snapshot_steps: Sequence[int] = (),
    perfusion_w_per_m3_k: float = 0.0,
    source_schedule: Optional[Callable[[float], float]] = None,
    t_init_expr: Optional[Callable] = None,
) -> TransientConductionResult:
    """The body of :func:`solve_transient_conduction`; ``perfusion_w_per_m3_k``
    of any sign (a negative value is `MAT-13`'s flipped-sign control only)."""
    for name, val in (("k", k), ("rho_c", rho_c)):
        if not _is_scalar(val):
            _check_positive_coefficient(mesh, name, val)
    if (
        (_is_scalar(k) and k <= 0.0)
        or h <= 0.0
        or (_is_scalar(rho_c) and rho_c <= 0.0)
        or dt <= 0.0
    ):
        raise ValueError("k, h, rho_c and dt must be positive")
    if not 0.0 <= theta <= 1.0:
        raise ValueError(f"theta must be in [0, 1], got {theta!r}")
    if startup_half_steps < 0 or startup_half_steps % 2:
        raise ValueError("startup_half_steps must be a non-negative even integer")
    if startup_half_steps // 2 > n_steps:
        raise ValueError("startup longer than the run")

    comm = mesh.comm
    V = fem.functionspace(mesh, ("Lagrange", degree))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)

    # A scalar k / rho_c becomes the same fem.Constant as before `MAT-18`
    # step 2(b); a UFL expression / fem.Function (DG0) is used as-is, and the
    # one rc_c object enters a_form, L_form and storage_form alike.
    k_c = _as_source(mesh, k)
    rc_c = _as_source(mesh, rho_c)
    q_c = _as_source(mesh, q)
    h_c = fem.Constant(mesh, default_scalar_type(h))
    tinf_c = fem.Constant(mesh, default_scalar_type(t_inf))
    dt_c = fem.Constant(mesh, default_scalar_type(dt))
    th_c = fem.Constant(mesh, default_scalar_type(theta))

    T_old = fem.Function(V, name="T_old")
    T_new = fem.Function(V, name="T")
    if t_init_expr is None:
        T_old.x.array[:] = t_init
        T_new.x.array[:] = t_init
    else:
        T_new.interpolate(t_init_expr)
        T_new.x.scatter_forward()
        T_old.x.array[:] = T_new.x.array

    perfused = perfusion_w_per_m3_k != 0.0
    w_c = fem.Constant(mesh, default_scalar_type(perfusion_w_per_m3_k)) if perfused else None

    def kh(a, z):
        out = k_c * ufl.inner(ufl.grad(a), ufl.grad(z)) * ufl.dx + h_c * ufl.inner(a, z) * ds(robin_tag)
        if perfused:
            out = out + w_c * ufl.inner(a, z) * ufl.dx
        return out

    a_form = fem.form(rc_c / dt_c * ufl.inner(u, v) * ufl.dx + th_c * kh(u, v))
    scheduled = source_schedule is not None
    s_c = fem.Constant(mesh, default_scalar_type(1.0)) if scheduled else None
    if scheduled:
        L_expr = (
            rc_c / dt_c * ufl.inner(T_old, v) * ufl.dx
            - (1.0 - th_c) * kh(T_old, v)
            + s_c * ufl.inner(q_c, v) * ufl.dx
            + ufl.inner(h_c * tinf_c, v) * ds(robin_tag)
        )
    else:
        L_expr = (
            rc_c / dt_c * ufl.inner(T_old, v) * ufl.dx
            - (1.0 - th_c) * kh(T_old, v)
            + ufl.inner(q_c, v) * ufl.dx
            + ufl.inner(h_c * tinf_c, v) * ds(robin_tag)
        )
    if perfused:
        L_expr = L_expr + ufl.inner(w_c * tinf_c, v) * ufl.dx
    L_form = fem.form(L_expr)
    storage_form = fem.form(rc_c * (T_new - T_old) / dt_c * ufl.dx)
    conv_form = fem.form(
        h_c * (th_c * T_new + (1.0 - th_c) * T_old - tinf_c) * ds(robin_tag)
    )
    perf_form = (
        fem.form(w_c * (th_c * T_new + (1.0 - th_c) * T_old - tinf_c) * ufl.dx)
        if perfused
        else None
    )
    src = float(np.real(comm.allreduce(
        fem.assemble_scalar(fem.form(q_c * ufl.dx(domain=mesh))), op=MPI.SUM
    )))
    # `MAT-21` step 19: with no source the old normalisation divides by zero;
    # normalise by the initial stored energy (per step, below) instead.
    e0 = None
    if src == 0.0:
        e0 = float(np.real(comm.allreduce(
            fem.assemble_scalar(fem.form(rc_c * abs(T_new - tinf_c) * ufl.dx)), op=MPI.SUM
        )))

    b = T_new.x.petsc_vec.duplicate()
    result = TransientConductionResult(temperature=T_new, source_power=src)
    snaps = set(int(s) for s in snapshot_steps)

    def _snap(step):
        if step in snaps:
            f = fem.Function(V, name=f"T_{step}")
            f.x.array[:] = T_new.x.array
            result.snapshots[step] = f

    _snap(0)

    def _run(dt_val, theta_val, n, step_of, t_mid_of=None):
        """Factor once for (dt_val, theta_val), then n back-substitutions."""
        dt_c.value = dt_val
        th_c.value = theta_val
        A = fem_petsc.assemble_matrix(a_form)
        A.assemble()
        prefix = f"fem_em_thermal_transient_{next(_PREFIX_COUNTER)}_"
        ksp = PETSc.KSP().create(comm)
        ksp.setOptionsPrefix(prefix)
        opts = PETSc.Options()
        opts[prefix + "ksp_type"] = "preonly"
        opts[prefix + "pc_type"] = "lu"
        opts[prefix + "pc_factor_mat_solver_type"] = "mumps"
        ksp.setFromOptions()
        ksp.setOperators(A)
        result.n_factorizations += 1
        try:
            for i in range(n):
                T_old.x.array[:] = T_new.x.array
                T_old.x.scatter_forward()
                if scheduled:
                    s_n = float(source_schedule(t_mid_of(i)))
                    s_c.value = s_n
                    result.source_multipliers.append(s_n)
                with b.localForm() as b_local:
                    b_local.set(0.0)
                fem_petsc.assemble_vector(b, L_form)
                b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
                ksp.solve(b, T_new.x.petsc_vec)
                T_new.x.scatter_forward()
                stor = comm.allreduce(fem.assemble_scalar(storage_form), op=MPI.SUM)
                conv = comm.allreduce(fem.assemble_scalar(conv_form), op=MPI.SUM)
                perf = 0.0
                if perfused:
                    perf = comm.allreduce(fem.assemble_scalar(perf_form), op=MPI.SUM)
                    result.perfused_power = float(np.real(perf))
                if e0 is not None:
                    res = abs(float(np.real(stor + conv + perf))) * dt_val / e0
                elif scheduled:
                    res = abs(float(np.real(stor + conv + perf)) - s_n * src) / abs(src)
                else:
                    res = abs(float(np.real(stor + conv + perf)) - src) / abs(src)
                result.balance_residuals.append(res)
                s = step_of(i)
                if s is not None:
                    _snap(s)
        finally:
            ksp.destroy()
            A.destroy()
            for key in ("ksp_type", "pc_type", "pc_factor_mat_solver_type"):
                opts.delValue(prefix + key)

    m = startup_half_steps
    if m:
        _run(0.5 * dt, 1.0, m, lambda i: (i + 1) // 2 if (i % 2) else None,
             lambda i: (i + 0.5) * (0.5 * dt))
    n_rest = n_steps - m // 2
    if n_rest:
        _run(dt, theta, n_rest, lambda i: m // 2 + i + 1,
             lambda i: (m // 2 + i + 0.5) * dt)

    b.destroy()
    result.max_balance_residual = max(result.balance_residuals, default=0.0)
    return result
