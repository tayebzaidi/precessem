"""`TH-5` step 3: a polynomial-graded PEC-backed PML on step 1a's TEM line.

Anchor: Jin 3e §9.6.1 (9.133) (metal-backed PML, |R| = exp(−2k₀s″L cosθ))
with s″L integrated over the §9.6.4 (9.165) profile s″ = s_max ξ^m:
Γ_c = −e^{−2jk₀d} exp(−2k₀ s_max d/(m + 1)) (a derivation, not a numbered
equation). Fixture: step 1a's TEM line (degree 1, vacuum, λ₀ = 1 m) extended
by a PML slab 0.9 < x < 2.4 m, PEC at x = 2.4 m; read by 1a's volume
projection over PROJ_N.

Review ruling 2026-10-22 03:00 (provisional; the 10-24 weekly ratifies or
reverts): the 1e-4 bar is asserted at 1a's k_fit read (the fixture's
dispersion leak removed); the k₀ read is printed beside it, never asserted
for the PML (its floor at λ₀/40 is ≈ 1.2e-4 for a reflectionless wave).
Normal incidence only.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem

from fem_em_solver.core.boundary_terms import (
    pml_reflection_pec_backed,
    pml_s_max_for_reflection,
    pml_slab_term,
)
from fem_em_solver.utils.constants import C_0, MU_0
from tests.validation.test_port_reaction_impedance import IDENTITY_TOLERANCE
from tests.validation.test_th5_abc_reflection import (
    K0,
    L_N,
    LAMBDA0,
    OMEGA,
    PROJ_N,
    TAG_DRIVE,
    TAG_TERM,
    TAG_Z,
    _assemble_complex,
    _box,
    _Projection,
    _tem,
    complex_only,
)

# The 2026-10-21 weekly: "≤ 1e-4 at λ₀/40 where the first-order ABC reads
# 9.5e-4 (step 1a's record)" (PROJECT_PLAN.md:4829).
PML_REFLECTION_BAR = 1e-4
PML_D = 1.5  # λ₀; 60 cells at λ₀/40; e^{−2jk₀d} = 1
PML_M = 2
R_ABSORB = 1e-6  # leg A
R_CLOSED = 1e-3  # leg C
S_ABSORB = pml_s_max_for_reflection(R_ABSORB, PML_D, PML_M, K0)
S_CLOSED = pml_s_max_for_reflection(R_CLOSED, PML_D, PML_M, K0)
ABC_K0_RECORD_40 = 9.499917e-4  # 20260930T094044Z_TH-5-step1a2-W1.log.gz:107
ETA0 = MU_0 * C_0


def _p(msg, comm=MPI.COMM_WORLD):
    if comm.rank == 0:
        print(msg, flush=True)


def test_pml_closed_form_numpy():
    """(Z) identities of the graded closed form, each ≤ IDENTITY_TOLERANCE."""
    worst = 0.0
    for r, s in ((R_ABSORB, S_ABSORB), (R_CLOSED, S_CLOSED)):
        g = pml_reflection_pec_backed(s, PML_D, PML_M, K0)
        worst = max(worst, abs(abs(g) / r - 1.0))
        # inversion
        s_back = pml_s_max_for_reflection(abs(g), PML_D, PML_M, K0)
        worst = max(worst, abs(s_back / s - 1.0))
        # σ form: s″ = σ/(ωε₀) ⇒ |R| = exp(−2η₀σd/(m+1))
        sigma = s * OMEGA / (MU_0 * C_0 * C_0)  # ε₀ = 1/(μ₀c₀²)
        r_sigma = np.exp(-2.0 * ETA0 * sigma * PML_D / (PML_M + 1))
        worst = max(worst, abs(r_sigma / r - 1.0))
        # angle: |Γ(θ)| = |Γ(0)|^{cosθ}
        for th in (0.3, 0.7, 1.1):
            gt = pml_reflection_pec_backed(s, PML_D, PML_M, K0, cos_theta=np.cos(th))
            worst = max(worst, abs(abs(gt) / abs(g) ** np.cos(th) - 1.0))
        # d = 1.5 λ₀ ⇒ Γ_c = −R
        worst = max(worst, abs(g + r) / r)
    g0 = pml_reflection_pec_backed(0.0, PML_D, PML_M, K0)
    worst = max(worst, abs(abs(g0) - 1.0))
    print(f"\n[TH-5 s3 (Z)] S_ABSORB = {S_ABSORB:.5f}, S_CLOSED = {S_CLOSED:.5f}, "
          f"worst identity = {worst:.3e}")
    assert abs(S_ABSORB - 2.19880) < 1e-5 and abs(S_CLOSED - 1.09940) < 1e-5
    assert worst <= IDENTITY_TOLERANCE, worst


def _solve_terms(msh, tags, degree, drive, pec_tags, extra):
    """1a's `_solve` (pec path) with extra bilinear terms added to `a`."""
    V = fem.functionspace(msh, ("N1curl", degree))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx - K0 * K0 * ufl.inner(u, v) * ufl.dx
    for t in extra:
        a = a + t(u, v)
    zero_vec = fem.Constant(msh, np.zeros(3, dtype=PETSc.ScalarType))
    L = ufl.inner(zero_vec, v) * ufl.dx

    fdim = msh.topology.dim - 1
    g = fem.Function(V)
    g.interpolate(drive)
    zero = fem.Function(V)
    zero.x.array[:] = 0.0
    bcs = [fem.dirichletbc(g, fem.locate_dofs_topological(V, fdim, tags.find(TAG_DRIVE)))]
    walls = list(pec_tags) + [TAG_TERM]
    wall_facets = np.unique(np.concatenate([tags.find(t) for t in walls]))
    bcs.append(fem.dirichletbc(zero, fem.locate_dofs_topological(V, fdim, wall_facets)))
    problem = LinearProblem(
        a, L, bcs=bcs,
        petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                       "pc_factor_mat_solver_type": "mumps"},
        petsc_options_prefix=f"fem_em_th5_pml_{degree}_",
    )
    e_h = problem.solve()
    e_h.x.scatter_forward()
    return e_h


def _pml_tem(comm, h, s_max, scan=True):
    t0 = time.perf_counter()
    length = L_N + PML_D
    w = 2.0 * h
    msh, tags, ncells = _box(comm, length, w, w, int(round(length / h)), 2, 2)
    drive = lambda x: np.vstack([np.zeros(x.shape[1]), np.zeros(x.shape[1]),
                                 np.ones(x.shape[1])]).astype(PETSc.ScalarType)
    term = pml_slab_term(msh, L_N, PML_D, s_max, PML_M, OMEGA)
    e_h = _solve_terms(msh, tags, 1, drive, [TAG_Z], [term])
    proj = _Projection(e_h[2], msh, PROJ_N[0], PROJ_N[1], K0)
    out = {"cells": ncells, "proj": proj.gamma(K0, L_N), "win_vol": proj.volume()}
    # (V) slab volume
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    chi.interpolate(lambda x: ((x[0] > L_N) & (x[0] < length)).astype(PETSc.ScalarType))
    out["slab_vol"] = _assemble_complex(fem.form(chi * ufl.dx), comm).real
    if scan:
        out["k_fit"], out["gamma_k_fit"], out["scan_interior"] = proj.k_fit(K0, L_N)
    out["seconds"] = time.perf_counter() - t0
    return out


@complex_only
@pytest.mark.integration
def test_pml_tem_line():
    comm = MPI.COMM_WORLD
    gc_a = pml_reflection_pec_backed(S_ABSORB, PML_D, PML_M, K0)
    gc_c = pml_reflection_pec_backed(S_CLOSED, PML_D, PML_M, K0)

    # Leg A ladder (printed), λ₀/40 asserted.
    ladder = {}
    for n in (20, 40, 80, 160):
        h = LAMBDA0 / n
        r = _pml_tem(comm, h, S_ABSORB)
        ladder[n] = r
        w = 2.0 * h
        rel_slab = abs(r["slab_vol"] / (PML_D * w * w) - 1.0)
        rel_win = abs(r["win_vol"] / (0.60 * w * w) - 1.0)
        _p(f"\n[TH-5 s3 A] lambda0/{n}: cells {r['cells']}, {r['seconds']:.1f} s; "
           f"k0 read |G-Gc| = {abs(r['proj'] - gc_a):.4e}; "
           f"k_fit-k0 = {r['k_fit'] - K0:+.4e} rad/m (interior {r['scan_interior']}), "
           f"k_fit read |G| = {abs(r['gamma_k_fit']):.4e}, |G-Gc| = {abs(r['gamma_k_fit'] - gc_a):.4e}; "
           f"slab vol rel {rel_slab:.2e}, window vol rel {rel_win:.2e}")
        assert r["cells"] == 24 * int(round((L_N + PML_D) / h)), r["cells"]
        assert rel_slab <= IDENTITY_TOLERANCE, rel_slab
        assert rel_win <= IDENTITY_TOLERANCE, rel_win
    ns = np.array([20, 40, 80])
    errs = np.array([abs(ladder[n]["gamma_k_fit"] - gc_a) for n in ns])
    rate = np.polyfit(np.log(1.0 / ns), np.log(errs), 1)[0]
    _p(f"[TH-5 s3 A] LSQ rate /20-/80 (k_fit read) = {rate:.3f}")
    a40 = ladder[40]
    assert a40["scan_interior"], a40["k_fit"] - K0
    err_a = abs(a40["gamma_k_fit"] - gc_a)

    # Leg C, λ₀/40.
    c40 = _pml_tem(comm, LAMBDA0 / 40, S_CLOSED)
    err_c = abs(c40["gamma_k_fit"] - gc_c)
    g_ungraded = -np.exp(-2.0 * K0 * S_CLOSED * PML_D)
    g_m3 = pml_reflection_pec_backed(S_CLOSED, PML_D, 3, K0)
    _p(f"[TH-5 s3 C] lambda0/40: Gc = {gc_c:.6e}; k_fit read G = {c40['gamma_k_fit']:.6e}, "
       f"|G| = {abs(c40['gamma_k_fit']):.6e}, |G-Gc| = {err_c:.4e} "
       f"(interior {c40['scan_interior']}); k0 read |G-Gc| = {abs(c40['proj'] - gc_c):.4e}")
    _p(f"[TH-5 s3 NC-grade] miss vs ungraded (9.133) {abs(c40['gamma_k_fit'] - g_ungraded):.4e}; "
       f"vs m=3 formula {abs(c40['gamma_k_fit'] - g_m3):.4e}")
    assert c40["scan_interior"], c40["k_fit"] - K0

    # NC-0 and NC-conj.
    z40 = _pml_tem(comm, LAMBDA0 / 40, 0.0, scan=False)
    j40 = _pml_tem(comm, LAMBDA0 / 40, -S_CLOSED, scan=False)
    _p(f"[TH-5 s3 NC-0] s_max = 0: k0 read G = {z40['proj']:.6e}, |G| = {abs(z40['proj']):.6f}")
    _p(f"[TH-5 s3 NC-conj] s_max = {-S_CLOSED:.5f}: k0 read |G| = {abs(j40['proj']):.4e}")

    # NC-ABC on the same line and read.
    abc = _tem(comm, LAMBDA0 / 40, 1, "abc", scan=True)
    ratio = abs(abc["gamma_k_fit"]) / abs(a40["gamma_k_fit"])
    _p(f"[TH-5 s3 NC-ABC] lambda0/40: k0 read |G| = {abs(abc['proj']):.6e} "
       f"(record {ABC_K0_RECORD_40:.6e}); k_fit-k0 = {abc['k_fit'] - K0:+.4e} "
       f"(interior {abc['scan_interior']}), k_fit read |G| = {abs(abc['gamma_k_fit']):.4e}; "
       f"|G_ABC|/|G_PML| at k_fit = {ratio:.1f}")
    _p(f"[TH-5 s3 summary] (A) |G-Gc| = {err_a:.4e}, (C) |G-Gc| = {err_c:.4e}, "
       f"bar {PML_REFLECTION_BAR:.1e}")

    assert abc["scan_interior"], abc["k_fit"] - K0
    assert abs(abc["proj"]) > PML_REFLECTION_BAR, abc["proj"]
    assert abs(abc["gamma_k_fit"]) > PML_REFLECTION_BAR, abc["gamma_k_fit"]
    assert abs(z40["proj"]) > PML_REFLECTION_BAR, z40["proj"]
    assert err_a <= PML_REFLECTION_BAR, err_a
    assert err_c <= PML_REFLECTION_BAR, err_c
