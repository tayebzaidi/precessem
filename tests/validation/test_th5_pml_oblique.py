"""`TH-5` step 3b: step 3's graded PEC-backed PML at oblique incidence.

Anchor: Jin 3e §9.6.1 (9.133) (metal-backed PML, |R| = exp(−2k₀s″L cosθ))
integrated over the §9.6.4 (9.165) profile s″ = s_max ξ^m (a derivation, not
a numbered equation; print-copy check is the operator's):
Γ_c(θ) = −e^{−2jk₀cosθ d} exp(−2k₀ cosθ s_max d/(m + 1)), so
|Γ_c(θ)| = R(0)^{cosθ}. At θ = 60°, d = 1.5 λ₀, R(0) = 1e-6: Γ_c = +1e-3.

Fixture: step 1a's TE10 guide at 60° (a = λ₀/(2 sin 60°), PEC y and z faces,
drive ẑ sin(πy/a) at x = 0, physical length L_O = 1.5 m) extended by a PML
slab L_O < x < L_O + PML_D, PEC at the far end; step 3's `pml_slab_term`
(m = 2, s_max = S_ABSORB from the normal-incidence design) and `_solve_terms`.
Read: 1a's `_Projection` with weight sin(πy/a) over FIT_O, asserted at the
k_fit read (β scanned ±5e-3 relative about k₀cosθ, 21 points; the read step 3
asserted at, ratified 10-24); the k₀cosθ read printed beside it.

§9 item 628 execution readings (2026-10-25 03:00 review): gate at degree 2,
degree 1 printed. Controls: (Cabc) first-order ABC |Γ| ≈ 1/3 (Jin (9.48)),
|Γ_ABC|/|Γ_PML| > 100; (NC-0) s_max = 0 reads |Γ| > 0.9. One angle, one
thickness, one profile.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core.boundary_terms import pml_reflection_pec_backed, pml_slab_term
from tests.validation.test_port_reaction_impedance import IDENTITY_TOLERANCE
from tests.validation.test_th5_abc_reflection import (
    A_O,
    BETA_O,
    FIT_O,
    K0,
    L_O,
    LAMBDA0,
    OMEGA,
    TAG_Y,
    TAG_Z,
    THETA_O,
    _box,
    _Projection,
    _solve,
    complex_only,
)
from tests.validation.test_th5_pml_reflection import (
    PML_D,
    PML_M,
    PML_REFLECTION_BAR,
    R_ABSORB,
    S_ABSORB,
    _solve_terms,
)

COS_O = float(np.cos(THETA_O))
N_GATE = 40  # λ₀/40
CABC_RATIO_BAR = 100.0  # ceiling (1/3)/1e-3 = 333
NC0_BAR = 0.9


def _p(msg, comm=MPI.COMM_WORLD):
    if comm.rank == 0:
        print(msg, flush=True)


def _drive(x):
    return np.vstack([np.zeros(x.shape[1]), np.zeros(x.shape[1]),
                      np.sin(np.pi * x[1] / A_O)]).astype(PETSc.ScalarType)


def _guide(comm, n_per_lambda, degree, s_max=None):
    """s_max given: PML-extended guide, PEC-backed. None: 1a's ABC guide."""
    t0 = time.perf_counter()
    h = LAMBDA0 / n_per_lambda
    ny = int(A_O / h)  # 1a's _te10
    wz = 2.0 * h
    if s_max is None:
        msh, tags, ncells = _box(comm, L_O, A_O, wz, int(round(L_O / h)), ny, 2)
        e_h = _solve(msh, tags, degree, _drive, [TAG_Y, TAG_Z], "abc")
    else:
        length = L_O + PML_D
        msh, tags, ncells = _box(comm, length, A_O, wz, int(round(length / h)), ny, 2)
        term = pml_slab_term(msh, L_O, PML_D, s_max, PML_M, OMEGA)
        e_h = _solve_terms(msh, tags, degree, _drive, [TAG_Y, TAG_Z], [term])
    x = ufl.SpatialCoordinate(msh)
    proj = _Projection(e_h[2], msh, FIT_O[0], FIT_O[1], BETA_O,
                       w=ufl.sin(ufl.pi * x[1] / A_O))
    out = {"cells": ncells, "proj": proj.gamma(BETA_O, L_O)}
    out["k_fit"], out["gamma_k_fit"], out["scan_interior"] = proj.k_fit(BETA_O, L_O)
    out["seconds"] = time.perf_counter() - t0
    return out


@complex_only
@pytest.mark.integration
def test_pml_oblique_te10():
    comm = MPI.COMM_WORLD
    gc = pml_reflection_pec_backed(S_ABSORB, PML_D, PML_M, K0, cos_theta=COS_O)
    gc0 = pml_reflection_pec_backed(S_ABSORB, PML_D, PML_M, K0)
    ident = abs(abs(gc) / abs(gc0) ** COS_O - 1.0)
    _p(f"\n[TH-5 s3b Gc] theta = 60 deg, cos = {COS_O:.15f}; Gc(theta) = {gc:.6e}, "
       f"|Gc| = {abs(gc):.6e} (R(0)^cos = {R_ABSORB ** COS_O:.6e}); "
       f"| |Gc|/|Gc0|^cos - 1 | = {ident:.3e}")

    r2 = _guide(comm, N_GATE, 2, S_ABSORB)
    r1 = _guide(comm, N_GATE, 1, S_ABSORB)
    for deg, r in ((2, r2), (1, r1)):
        _p(f"[TH-5 s3b B] degree {deg}, lambda0/{N_GATE}: cells {r['cells']}, "
           f"{r['seconds']:.1f} s; k_fit-beta = {r['k_fit'] - BETA_O:+.4e} rad/m "
           f"(interior {r['scan_interior']}); k_fit read G = {r['gamma_k_fit']:.6e}, "
           f"|G| = {abs(r['gamma_k_fit']):.6e}, |G-Gc| = {abs(r['gamma_k_fit'] - gc):.4e}; "
           f"k0cos read G = {r['proj']:.6e}, |G-Gc| = {abs(r['proj'] - gc):.4e}")

    abc = _guide(comm, N_GATE, 2, None)
    nc0 = _guide(comm, N_GATE, 2, 0.0)
    ratio = abs(abc["gamma_k_fit"]) / abs(r2["gamma_k_fit"])
    _p(f"[TH-5 s3b Cabc] degree 2 ABC: {abc['seconds']:.1f} s; k_fit read G = "
       f"{abc['gamma_k_fit']:.6e} (|G| {abs(abc['gamma_k_fit']):.6e}, interior "
       f"{abc['scan_interior']}); k0cos read G = {abc['proj']:.6e}; "
       f"|G_ABC|/|G_PML| (k_fit) = {ratio:.1f}")
    _p(f"[TH-5 s3b NC-0] s_max = 0: k_fit read |G| = {abs(nc0['gamma_k_fit']):.6f}, "
       f"k0cos read |G| = {abs(nc0['proj']):.6f}")
    err = abs(r2["gamma_k_fit"] - gc)
    _p(f"[TH-5 s3b summary] (B) degree 2 |G-Gc| = {err:.4e}, bar {PML_REFLECTION_BAR:.1e}; "
       f"degree 1 {abs(r1['gamma_k_fit'] - gc):.4e} (printed)")

    assert abs(gc - 1e-3) <= 1e-3 * 1e-9, gc
    assert ident <= IDENTITY_TOLERANCE, ident
    assert r2["scan_interior"], r2["k_fit"] - BETA_O
    assert ratio > CABC_RATIO_BAR, ratio
    assert abs(nc0["gamma_k_fit"]) > NC0_BAR, nc0["gamma_k_fit"]
    assert err <= PML_REFLECTION_BAR, err
