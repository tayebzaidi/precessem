"""`TH-5` step 3c: the graded PML's angular law at 45°, on the closed leg.

Anchor: Jin 3e §9.6.1 (9.133) (metal-backed PML, |R| = exp(−2k₀s″L cosθ))
integrated over the §9.6.4 (9.165) profile s″ = s_max ξ^m (a derivation, not
a numbered equation; print-copy check is the operator's):
Γ_c(θ) = −e^{−2jk₀cosθ d} exp(−2k₀ cosθ s_max d/(m + 1)), so
|Γ_c(θ)| = R(0)^{cosθ}. Closed leg R(0) = R_CLOSED = 1e-3 (step 3's leg C):
|Γ_c(45°)| = 1e-3^{0.7071} ≈ 7.6e-3, 76× PML_REFLECTION_BAR, phase
−2k₀ cos45° d (d = 1.5 λ₀) not a multiple of π, so Γ_c is complex.

Fixture: step 1a's TE10 guide geometry rebuilt at θ = 45° (a = λ₀/(2 sin θ),
β = k₀ cos θ, computed locally — never the 60° module constants), extended by
step 3's `pml_slab_term` slab (m = 2, d = 1.5 λ₀, s_max = S_CLOSED), PEC at
the far end; step 3's `_solve_terms`. Read: 1a's `_Projection` weighted by
sin(πy/a) over FIT_O, asserted at the k_fit read (β scanned ±5e-3 relative,
21 points), the k₀cosθ read printed beside it.

§9 item 644 (2026-10-25 18:00 review). Asserted: (Gc45) the magnitude
identity; (B45) degree 2, λ₀/40, |Γ_h − Γ_c(45°)| ≤ PML_REFLECTION_BAR;
(Law) the normal-incidence law misses, |Γ_h − Γ_c(cos θ = 1)| > 10 × bar
(ceiling ≥ 66×). Printed only: degree 1, the 60° law on this leg, and (Cabc)
(ceiling ≈ 22.6 on this leg, under CABC_RATIO_BAR by construction).
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
    FIT_O,
    K0,
    L_O,
    LAMBDA0,
    OMEGA,
    TAG_Y,
    TAG_Z,
    _box,
    _Projection,
    _solve,
    complex_only,
)
from tests.validation.test_th5_pml_oblique import CABC_RATIO_BAR, N_GATE
from tests.validation.test_th5_pml_reflection import (
    PML_D,
    PML_M,
    PML_REFLECTION_BAR,
    R_CLOSED,
    S_CLOSED,
    _solve_terms,
)

THETA_45 = np.deg2rad(45.0)
COS_45 = float(np.cos(THETA_45))
LAW_FACTOR = 10.0  # (Law): miss the normal-incidence law by > 10 × bar (ceiling ≥ 66×)


def _p(msg, comm=MPI.COMM_WORLD):
    if comm.rank == 0:
        print(msg, flush=True)


def _closed_forms():
    gc = pml_reflection_pec_backed(S_CLOSED, PML_D, PML_M, K0, cos_theta=COS_45)
    gc0 = pml_reflection_pec_backed(S_CLOSED, PML_D, PML_M, K0)
    gc60 = pml_reflection_pec_backed(S_CLOSED, PML_D, PML_M, K0, cos_theta=0.5)
    return gc, gc0, gc60


def test_pml_angle_closed_forms_numpy():
    """(Gc45) identity and the (Law) ceiling, no mesh."""
    gc, gc0, _ = _closed_forms()
    ident = abs(abs(gc) / abs(gc0) ** COS_45 - 1.0)
    target = R_CLOSED ** COS_45
    _p(f"\n[TH-5 s3c W0] Gc(45) = {gc:.6e}, |Gc| = {abs(gc):.6e} (R^cos = {target:.6e}); "
       f"ident {ident:.3e}; Gc(0) = {gc0:.6e}; |Gc45 - Gc0| = {abs(gc - gc0):.4e} "
       f"({abs(gc - gc0) / PML_REFLECTION_BAR:.1f}x bar)")
    assert ident <= IDENTITY_TOLERANCE, ident
    assert abs(abs(gc) - target) <= 1e-12 * target, (abs(gc), target)
    assert abs(gc.imag) > 1e-2 * abs(gc), gc  # complex: phase not a multiple of π
    assert abs(gc - gc0) > LAW_FACTOR * PML_REFLECTION_BAR, abs(gc - gc0)


def _guide_at(comm, theta, n_per_lambda, degree, s_max):
    """s_max given: PML-extended 45° guide, PEC-backed. None: ABC guide."""
    t0 = time.perf_counter()
    a = LAMBDA0 / (2.0 * np.sin(theta))
    beta = K0 * np.cos(theta)
    h = LAMBDA0 / n_per_lambda
    ny = int(a / h)
    wz = 2.0 * h

    def drive(x):
        return np.vstack([np.zeros(x.shape[1]), np.zeros(x.shape[1]),
                          np.sin(np.pi * x[1] / a)]).astype(PETSc.ScalarType)

    if s_max is None:
        msh, tags, ncells = _box(comm, L_O, a, wz, int(round(L_O / h)), ny, 2)
        e_h = _solve(msh, tags, degree, drive, [TAG_Y, TAG_Z], "abc")
    else:
        length = L_O + PML_D
        msh, tags, ncells = _box(comm, length, a, wz, int(round(length / h)), ny, 2)
        term = pml_slab_term(msh, L_O, PML_D, s_max, PML_M, OMEGA)
        e_h = _solve_terms(msh, tags, degree, drive, [TAG_Y, TAG_Z], [term])
    x = ufl.SpatialCoordinate(msh)
    proj = _Projection(e_h[2], msh, FIT_O[0], FIT_O[1], beta,
                       w=ufl.sin(ufl.pi * x[1] / a))
    out = {"cells": ncells, "ny": ny, "beta": beta, "proj": proj.gamma(beta, L_O)}
    out["k_fit"], out["gamma_k_fit"], out["scan_interior"] = proj.k_fit(beta, L_O)
    out["seconds"] = time.perf_counter() - t0
    return out


@complex_only
@pytest.mark.integration
def test_pml_angle_45_closed_leg():
    comm = MPI.COMM_WORLD
    gc, gc0, gc60 = _closed_forms()
    ident = abs(abs(gc) / abs(gc0) ** COS_45 - 1.0)
    _p(f"\n[TH-5 s3c Gc45] theta = 45 deg, cos = {COS_45:.15f}; Gc(45) = {gc:.6e}, "
       f"|Gc| = {abs(gc):.6e} (R_CLOSED^cos = {R_CLOSED ** COS_45:.6e}); "
       f"| |Gc|/|Gc0|^cos - 1 | = {ident:.3e}; Gc(0) = {gc0:.6e}; Gc(60) = {gc60:.6e}")

    r2 = _guide_at(comm, THETA_45, N_GATE, 2, S_CLOSED)
    r1 = _guide_at(comm, THETA_45, N_GATE, 1, S_CLOSED)
    for deg, r in ((2, r2), (1, r1)):
        _p(f"[TH-5 s3c B45] degree {deg}, lambda0/{N_GATE}: cells {r['cells']} (ny {r['ny']}), "
           f"{r['seconds']:.1f} s; k_fit-beta = {r['k_fit'] - r['beta']:+.4e} rad/m "
           f"(interior {r['scan_interior']}); k_fit read G = {r['gamma_k_fit']:.6e}, "
           f"|G| = {abs(r['gamma_k_fit']):.6e}, |G-Gc45| = {abs(r['gamma_k_fit'] - gc):.4e}; "
           f"k0cos read G = {r['proj']:.6e}, |G-Gc45| = {abs(r['proj'] - gc):.4e}")

    g = r2["gamma_k_fit"]
    err = abs(g - gc)
    law0 = abs(g - gc0)
    law60 = abs(g - gc60)
    _p(f"[TH-5 s3c Law] degree 2 k_fit read: |G-Gc(cos=1)| = {law0:.4e} "
       f"({law0 / PML_REFLECTION_BAR:.1f}x bar; asserted > {LAW_FACTOR:.0f}x); "
       f"|G-Gc(60)| = {law60:.4e} (printed)")

    abc = _guide_at(comm, THETA_45, N_GATE, 2, None)
    ratio = abs(abc["gamma_k_fit"]) / abs(g)
    g_abc_c = (COS_45 - 1.0) / (COS_45 + 1.0)
    _p(f"[TH-5 s3c Cabc] degree 2 ABC (printed only): {abc['seconds']:.1f} s; k_fit read G = "
       f"{abc['gamma_k_fit']:.6e} (interior {abc['scan_interior']}), Jin (9.48) "
       f"{g_abc_c:.6e}; k0cos read G = {abc['proj']:.6e}; |G_ABC|/|G_PML| = {ratio:.2f} "
       f"(CABC_RATIO_BAR {CABC_RATIO_BAR:.0f} not applied; ceiling ~22.6)")
    _p(f"[TH-5 s3c summary] (B45) degree 2 |G-Gc45| = {err:.4e}, bar {PML_REFLECTION_BAR:.1e}; "
       f"degree 1 {abs(r1['gamma_k_fit'] - gc):.4e} (printed); (Law) {law0:.4e}")

    assert ident <= IDENTITY_TOLERANCE, ident
    assert abs(abs(gc) - R_CLOSED ** COS_45) <= 1e-12 * R_CLOSED ** COS_45, abs(gc)
    assert r2["scan_interior"], r2["k_fit"] - r2["beta"]
    assert law0 > LAW_FACTOR * PML_REFLECTION_BAR, law0
    assert err <= PML_REFLECTION_BAR, err
