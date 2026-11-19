"""`TH-5` step 1a: the first-order ABC as a surface term, on Jin's reflection.

Anchor: Jin 3e §9.2 (9.44) at s = 1 / §9.3 (9.54); plane-wave TE reflection of
the first-order ABC R(θ) = (cosθ − 1)/(cosθ + 1) (9.48), the first-order Padé
(9.4) of §9.1.1 — 0 at θ = 0, −1/3 at 60°, |R| → 1 at grazing. Re-derived here
by the mode-impedance route Γ = (η₀ − Z_TE)/(η₀ + Z_TE), Z_TE = η₀/cosθ.

FEM fixtures (degree 1, vacuum, f = c₀ Hz so λ₀ = 1 m, e^{+jωt}):

(N) TEM line: PEC z faces, natural y faces, L = 0.9 m, W = 2h; drive ẑ at
    x = 0, ABC on x = L. Γ = (B/A)e^{2jk₀L}, with (A, B) the **volume
    least-squares projection** of E_z onto {e^{−jk₀x}, e^{+jk₀x}} over the
    whole cell columns of x ∈ [0.15, 0.75] m (attempt 2, §9 item 94). The
    41-point fit on [0.135, 0.765] m of attempt 1 is printed beside it: its
    points step 0.315h … 2.52h, each at its own in-cell offset, and read a
    flat ≈ 1.3–1.5e-3 floor the projection does not see.
(O) TE10 guide at 60°: PEC y and z faces, a = λ₀/(2 sin 60°), β = π rad/m,
    L = 1.5 m; drive ẑ sin(πy/a). (O1) Jin's exact condition γ = jk₀cosθ gives
    Γ = 0 (premise); (O2) the first-order ABC gives −1/3.

(P) Projection premises, no solve: a planted two-wave field is recovered to
    1e-10, and ∫χ dV = 0.60·(2h)² to 1e-12 at every (N) rung.

(O)'s gate is the 41-point fit, unchanged (its points sit on lattice planes);
its projection is printed.

Controls (PEC / PMC termination) are printed beside the measured value, never
asserted (rule (e): new fixture).
"""

from __future__ import annotations

import os

import numpy as np
import pytest
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core.boundary_terms import (
    abc_first_order_reflection_te,
    first_order_abc_term,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.constants import C_0

IS_COMPLEX = np.issubdtype(np.dtype(PETSc.ScalarType), np.complexfloating)
REQUIRE_COMPLEX = os.environ.get("FEM_EM_REQUIRE_COMPLEX", "") == "1"
complex_only = pytest.mark.skipif(
    not IS_COMPLEX and not REQUIRE_COMPLEX,
    reason=f"TH-5 needs the complex build; this one is {np.dtype(PETSc.ScalarType)}",
)

LAMBDA0 = 1.0
FREQUENCY_HZ = C_0 / LAMBDA0
OMEGA = 2.0 * np.pi * FREQUENCY_HZ
K0 = 2.0 * np.pi / LAMBDA0

TAG_DRIVE, TAG_TERM, TAG_Y, TAG_Z = 1, 2, 3, 4

# (N) TEM line
L_N = 0.9
FIT_N = (0.135, 0.765)
# Projection window: FIT_N snapped inward to the λ₀/20 planes (whole cell
# columns at every rung: 12 / 24 / 48 / 96).
PROJ_N = (0.15, 0.75)
PLANTED_GAMMA = 0.01 - 0.02j
# Attempt 1's 41-point (N) ladder (`20260930T022458Z_TH-5-step1a-W1.log:95–102`).
POINT_FIT_RECORD = (4.916e-3, 2.018e-3, 1.531e-3, 1.296e-3)
# (O) TE10 at 60 degrees
THETA_O = np.deg2rad(60.0)
A_O = LAMBDA0 / (2.0 * np.sin(THETA_O))
L_O = 1.5
BETA_O = K0 * np.cos(THETA_O)  # = pi rad/m
FIT_O = (0.25, 1.25)
O_BAND = 3e-3  # (O1) / (O2) band, hoisted by value for `th:18` (EX-143)


# ----------------------------------------------------------------- numpy ---


def test_numpy_reflection_anchors():
    """R(1) = 0, R(0.5) = −1/3, grazing → −1, mode-impedance route agrees."""
    assert abc_first_order_reflection_te(1.0) == 0.0
    assert abs(abc_first_order_reflection_te(0.5) - (-1.0 / 3.0)) <= 1e-15
    r_near_grazing = abc_first_order_reflection_te(np.cos(np.deg2rad(89.9)))
    assert -1.0 < r_near_grazing < 0.0
    assert abs(abc_first_order_reflection_te(1e-12) - (-1.0)) < 1e-11

    eta0 = 376.730313461
    theta = np.linspace(0.0, np.deg2rad(89.0), 50)
    z_te = eta0 / np.cos(theta)
    gamma_mode = (eta0 - z_te) / (eta0 + z_te)
    worst = float(np.max(np.abs(gamma_mode - abc_first_order_reflection_te(np.cos(theta)))))
    print(
        f"\n[TH-5 1a numpy] R(0.5) = {abc_first_order_reflection_te(0.5):.17f}; "
        f"R(cos 89.9deg) = {r_near_grazing:.6f}; mode-impedance max diff = {worst:.3e}"
    )
    assert worst <= 1e-15


# ------------------------------------------------------------------- FEM ---


def _box(comm, length, wy, wz, nx, ny, nz):
    msh = dmesh.create_box(
        comm,
        [np.zeros(3), np.array([length, wy, wz])],
        [nx, ny, nz],
        cell_type=dmesh.CellType.tetrahedron,
    )
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    tol = 1e-9
    marks = [
        (TAG_DRIVE, lambda x: np.isclose(x[0], 0.0, atol=tol)),
        (TAG_TERM, lambda x: np.isclose(x[0], length, atol=tol)),
        (TAG_Y, lambda x: np.isclose(x[1], 0.0, atol=tol) | np.isclose(x[1], wy, atol=tol)),
        (TAG_Z, lambda x: np.isclose(x[2], 0.0, atol=tol) | np.isclose(x[2], wz, atol=tol)),
    ]
    idx, val = [], []
    for tag, fn in marks:
        f = dmesh.locate_entities_boundary(msh, fdim, fn)
        idx.append(f)
        val.append(np.full(f.size, tag, dtype=np.int32))
    idx = np.concatenate(idx)
    val = np.concatenate(val)
    order = np.argsort(idx, kind="stable")
    idx, val = idx[order], val[order]
    # A facet on two planes (none here: faces are disjoint in facet sets) —
    # keep the first occurrence so meshtags indices stay unique.
    uniq, first = np.unique(idx, return_index=True)
    tags = dmesh.meshtags(msh, fdim, uniq, val[first])
    ncells = comm.allreduce(msh.topology.index_map(msh.topology.dim).size_local, op=MPI.SUM)
    return msh, tags, ncells


def _solve(msh, tags, degree, drive, pec_tags, termination, cos_theta=1.0):
    """termination: 'abc' | 'exact' | 'pec' | 'pmc'."""
    V = fem.functionspace(msh, ("N1curl", degree))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx - K0 * K0 * ufl.inner(u, v) * ufl.dx
    if termination == "abc":
        a = a + first_order_abc_term(msh, tags, TAG_TERM, OMEGA)(u, v)
    elif termination == "exact":
        a = a + first_order_abc_term(msh, tags, TAG_TERM, OMEGA, cos_theta=cos_theta)(u, v)
    zero_vec = fem.Constant(msh, np.zeros(3, dtype=PETSc.ScalarType))
    L = ufl.inner(zero_vec, v) * ufl.dx

    fdim = msh.topology.dim - 1
    g = fem.Function(V)
    g.interpolate(drive)
    zero = fem.Function(V)
    zero.x.array[:] = 0.0
    bcs = [fem.dirichletbc(g, fem.locate_dofs_topological(V, fdim, tags.find(TAG_DRIVE)))]
    walls = list(pec_tags) + ([TAG_TERM] if termination == "pec" else [])
    wall_facets = np.unique(np.concatenate([tags.find(t) for t in walls]))
    # Wall dofs shared with the drive face carry zero tangential drive (ẑ on
    # z-faces is normal; sin(πy/a) vanishes on y-faces), so the order is moot.
    bcs.append(fem.dirichletbc(zero, fem.locate_dofs_topological(V, fdim, wall_facets)))
    problem = LinearProblem(
        a,
        L,
        bcs=bcs,
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
        petsc_options_prefix=f"fem_em_th5_{termination}_{degree}_",
    )
    e_h = problem.solve()
    e_h.x.scatter_forward()
    return e_h


def _fit_gamma(e_h, comm, x_range, y0, z0, beta, length, scale=1.0, npts=41):
    xs = np.linspace(x_range[0], x_range[1], npts)
    pts = np.column_stack([xs, np.full_like(xs, y0), np.full_like(xs, z0)])
    values, valid = evaluate_vector_field_parallel(e_h, pts, comm)
    gamma = None
    if comm.rank == 0:
        assert np.all(valid), "fit points not all located"
        ez = np.asarray(values)[:, 2] / scale
        M = np.column_stack([np.exp(-1j * beta * xs), np.exp(1j * beta * xs)])
        (A, B), *_ = np.linalg.lstsq(M, ez, rcond=None)
        gamma = complex(B / A * np.exp(2j * beta * length))
    return comm.bcast(gamma, root=0)


def _assemble_complex(compiled, comm):
    """Rank-summed complex scalar (no ``np.real``: the projection is complex)."""
    return complex(comm.allreduce(fem.assemble_scalar(compiled), op=MPI.SUM))


class _Projection:
    """Volume LSQ projection of ``ez`` onto w·e^{∓jβx} over whole cell columns.

    χ is a DG0 indicator set in numpy at cell midpoints (never a UFL
    comparison).  M[p][q] = ∫χ φ_q conj(φ_p), b[p] = ∫χ ez conj(φ_p): M's
    transpose is its conjugate and reads |Γ| ≈ 0.2 on a pure forward wave.
    β lives in a Constant, so a scan over β reuses the compiled forms.
    """

    def __init__(self, ez, msh, x_lo, x_hi, beta, w=1.0):
        self.comm = msh.comm
        dg0 = fem.functionspace(msh, ("DG", 0))
        chi = fem.Function(dg0)
        chi.interpolate(lambda x: ((x[0] > x_lo) & (x[0] < x_hi)).astype(PETSc.ScalarType))
        self.kc = fem.Constant(msh, PETSc.ScalarType(beta))
        j = fem.Constant(msh, PETSc.ScalarType(1j))
        x = ufl.SpatialCoordinate(msh)
        c, s = ufl.cos(self.kc * x[0]), ufl.sin(self.kc * x[0])
        phi = (w * (c - j * s), w * (c + j * s))
        dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 6})
        self._m = [
            [fem.form(chi * ufl.inner(phi[q], phi[p]) * dx) for q in range(2)] for p in range(2)
        ]
        self._b = [fem.form(chi * ufl.inner(ez, phi[p]) * dx) for p in range(2)]
        self._e2 = fem.form(chi * ufl.inner(ez, ez) * dx)
        self._vol = fem.form(chi * dx)

    def volume(self):
        return _assemble_complex(self._vol, self.comm).real

    def solve(self, beta):
        """(A, B) and the squared residual ∫χ|ez − Aφ₀ − Bφ₁|² at ``beta``."""
        self.kc.value = PETSc.ScalarType(beta)
        M = np.array(
            [[_assemble_complex(self._m[p][q], self.comm) for q in range(2)] for p in range(2)]
        )
        b = np.array([_assemble_complex(self._b[p], self.comm) for p in range(2)])
        ab = np.linalg.solve(M, b)
        e2 = _assemble_complex(self._e2, self.comm).real
        return ab, float(e2 - np.real(np.vdot(b, ab)))

    def gamma(self, beta, length):
        (A, B), _ = self.solve(beta)
        return complex(B / A * np.exp(2j * beta * length))

    def k_fit(self, k0, length, half_width=5e-3, n=21):
        """Residual-minimising k: 21-point scan, parabolic vertex; Γ there."""
        ks = k0 * (1.0 + np.linspace(-half_width, half_width, n))
        res = np.array([self.solve(k)[1] for k in ks])
        i = int(np.clip(np.argmin(res), 1, n - 2))
        a2, a1, _ = np.polyfit(ks[i - 1 : i + 2] - ks[i], res[i - 1 : i + 2], 2)
        k = float(ks[i] - a1 / (2.0 * a2))
        return k, self.gamma(k, length), bool(0 < int(np.argmin(res)) < n - 1)


def _project_gamma(ez, msh, x_lo, x_hi, beta, length, w=1.0):
    return _Projection(ez, msh, x_lo, x_hi, beta, w).gamma(beta, length)


def _tem_box(comm, h):
    nx = int(round(L_N / h))
    w = 2.0 * h
    return _box(comm, L_N, w, w, nx, 2, 2) + (w,)


def _tem(comm, h, degree, termination, scan=False):
    """One solve, both readings: the 41-point fit and the volume projection."""
    msh, tags, ncells, w = _tem_box(comm, h)
    drive = lambda x: np.vstack(
        [np.zeros(x.shape[1]), np.zeros(x.shape[1]), np.ones(x.shape[1])]
    ).astype(PETSc.ScalarType)
    e_h = _solve(msh, tags, degree, drive, [TAG_Z], termination)
    out = {"cells": ncells, "fit": _fit_gamma(e_h, comm, FIT_N, 0.37 * w, 0.41 * w, K0, L_N)}
    proj = _Projection(e_h[2], msh, PROJ_N[0], PROJ_N[1], K0)
    out["proj"] = proj.gamma(K0, L_N)
    if scan:
        out["k_fit"], out["gamma_k_fit"], out["scan_interior"] = proj.k_fit(K0, L_N)
    return out


def _te10(comm, n_per_lambda, degree, termination, cos_theta=1.0, return_fields=False):
    h = LAMBDA0 / n_per_lambda
    nx = int(round(L_O / h))
    ny = int(A_O / h)
    wz = 2.0 * h
    msh, tags, ncells = _box(comm, L_O, A_O, wz, nx, ny, 2)
    drive = lambda x: np.vstack(
        [np.zeros(x.shape[1]), np.zeros(x.shape[1]), np.sin(np.pi * x[1] / A_O)]
    ).astype(PETSc.ScalarType)
    e_h = _solve(msh, tags, degree, drive, [TAG_Y, TAG_Z], termination, cos_theta)
    y0 = 0.437 * A_O
    gamma = _fit_gamma(e_h, comm, FIT_O, y0, 0.41 * wz, BETA_O, L_O, scale=np.sin(np.pi * y0 / A_O))
    x = ufl.SpatialCoordinate(msh)
    gamma_proj = _project_gamma(
        e_h[2], msh, FIT_O[0], FIT_O[1], BETA_O, L_O, w=ufl.sin(ufl.pi * x[1] / A_O)
    )
    if return_fields:
        return gamma, ncells, (nx, ny, 2), gamma_proj, e_h, msh, tags
    return gamma, ncells, (nx, ny, 2), gamma_proj


@complex_only
@pytest.mark.integration
def test_projection_premises():
    """(P), no solve: a planted two-wave field is recovered; ∫χ dV is whole columns."""
    comm = MPI.COMM_WORLD
    worst_vol = 0.0
    for n in (20, 40, 80, 160):
        h = LAMBDA0 / n
        msh, _tags, _nc, w = _tem_box(comm, h)
        x = ufl.SpatialCoordinate(msh)
        k = fem.Constant(msh, PETSc.ScalarType(K0))
        j = fem.Constant(msh, PETSc.ScalarType(1j))
        back = fem.Constant(msh, PETSc.ScalarType(PLANTED_GAMMA * np.exp(-2j * K0 * L_N)))
        c, s = ufl.cos(k * x[0]), ufl.sin(k * x[0])
        planted = (c - j * s) + back * (c + j * s)
        proj = _Projection(planted, msh, PROJ_N[0], PROJ_N[1], K0)
        vol = proj.volume()
        rel = abs(vol / (0.60 * w * w) - 1.0)
        worst_vol = max(worst_vol, rel)
        if comm.rank == 0:
            print(
                f"\n[TH-5 1a2 (P)] lambda0/{n}: int chi dV = {vol:.15e}, "
                f"|./(0.60 (2h)^2) - 1| = {rel:.3e}"
            )
        assert rel <= 1e-12, (n, rel)
        if n == 20:
            g = proj.gamma(K0, L_N)
            if comm.rank == 0:
                print(
                    f"[TH-5 1a2 (P)] planted Gamma = {PLANTED_GAMMA}, recovered {g:.12e}, "
                    f"|diff| = {abs(g - PLANTED_GAMMA):.3e}"
                )
            assert abs(g - PLANTED_GAMMA) <= 1e-10, g


@complex_only
@pytest.mark.integration
def test_tem_line_abc_ladder():
    """(N), by projection: |Γ| ≤ 1e-3 at λ₀/80 and λ₀/160, decreasing, rate ≥ 1.5."""
    comm = MPI.COMM_WORLD
    ns = [20, 40, 80, 160]
    hs = [LAMBDA0 / n for n in ns]
    gammas, points = [], []
    abc80 = None
    for n, h in zip(ns, hs):
        r = _tem(comm, h, 1, "abc", scan=(n == 80))
        if n == 80:
            abc80 = r
        gammas.append(abs(r["proj"]))
        points.append(abs(r["fit"]))
        if comm.rank == 0:
            print(
                f"\n[TH-5 1a2 (N)] h = lambda0/{n} ({r['cells']} cells): projection Gamma = "
                f"{r['proj']:.6e}, |Gamma| = {abs(r['proj']):.6e}; 41-point fit |Gamma| = "
                f"{abs(r['fit']):.6e} (attempt 1: {POINT_FIT_RECORD[ns.index(n)]:.3e})"
            )
    rate = float(np.polyfit(np.log(hs), np.log(gammas), 1)[0])
    rate_pts = float(np.polyfit(np.log(hs), np.log(points), 1)[0])
    d2 = _tem(comm, LAMBDA0 / 40, 2, "abc")
    # Controls (printed only, rule (e)).
    pec = _tem(comm, LAMBDA0 / 80, 1, "pec", scan=True)
    pmc = _tem(comm, LAMBDA0 / 80, 1, "pmc")
    if comm.rank == 0:
        print(
            f"[TH-5 1a2 (N)] LSQ rate of |Gamma| in h: projection {rate:.4f} "
            f"(predicted ~2), 41-point fit {rate_pts:.4f} (attempt 1: 0.617)"
        )
        print(
            "[TH-5 1a2 (N)] predicted projection ladder 2.3e-3 / 5.7e-4 / 1.4e-4 / 3.5e-5 "
            "(x/ 3); sampling hypothesis refuted if >= 9.1e-4 at lambda0/160, confirmed if "
            "<= 4.3e-4 / <= 1.1e-4 at lambda0/80 / /160"
        )
        print(
            f"[TH-5 1a2 (N)] degree 2 at lambda0/40 ({d2['cells']} cells): projection |Gamma| = "
            f"{abs(d2['proj']):.6e}, 41-point fit {abs(d2['fit']):.6e} (printed)"
        )
        for name, c in (("PEC", pec), ("PMC", pmc)):
            g0 = c["proj"] * np.exp(-2j * K0 * L_N)
            print(
                f"[TH-5 1a2 (N) control {name} @ lambda0/80] projection Gamma = {c['proj']:.6f}, "
                f"|1+Gamma(0)| = {abs(1 + g0):.4f}, "
                f"factor |Gamma_{name}|/|Gamma_ABC| = {abs(c['proj']) / gammas[2]:.1f}; "
                f"41-point fit Gamma = {c['fit']:.6f}"
            )
        print(
            f"[TH-5 1a2 (N) control PEC] Im Gamma by projection = {pec['proj'].imag:.4e} "
            f"(predicted 3.2e-3 x/ 3 if dispersion; point fit read 0.0272)"
        )
        for name, c in (("ABC", abc80), ("PEC", pec)):
            print(
                f"[TH-5 1a2 (N) k_fit {name} @ lambda0/80] k_fit - k0 = {c['k_fit'] - K0:+.4e} rad/m "
                f"(predicted +1.8e-3; vertex interior: {c['scan_interior']}), "
                f"|Gamma| at k_fit = {abs(c['gamma_k_fit']):.6e} (printed, never the gate)"
            )
    assert gammas[2] <= 1e-3 and gammas[3] <= 1e-3, gammas
    assert all(gammas[i + 1] < gammas[i] for i in range(3)), gammas
    assert rate >= 1.5, rate


@complex_only
@pytest.mark.integration
def test_te10_sixty_degrees():
    """(O1) exact-angle condition |Γ| ≤ 3e-3; (O2) ABC |Γ − (−1/3)| ≤ 3e-3."""
    comm = MPI.COMM_WORLD
    g_exact, nc, shape, _ = _te10(comm, 80, 1, "exact", cos_theta=np.cos(THETA_O))
    g_abc, _, _, p_abc = _te10(comm, 80, 1, "abc")
    g_abc40, _, _, p_abc40 = _te10(comm, 40, 1, "abc")
    g_abc40_d2, _, _, p_abc40_d2 = _te10(comm, 40, 2, "abc")
    g_pec, _, _, _ = _te10(comm, 80, 1, "pec")
    r_ref = float(abc_first_order_reflection_te(np.cos(THETA_O)))
    if comm.rank == 0:
        print(f"\n[TH-5 1a2 (O)] mesh {shape}, {nc} cells, a = {A_O:.5f} m, beta = {BETA_O:.6f}")
        print(f"[TH-5 1a2 (O1)] exact-angle Gamma = {g_exact:.6e}, |Gamma| = {abs(g_exact):.6e}")
        print(
            f"[TH-5 1a2 (O2)] ABC Gamma = {g_abc:.6e}, |Gamma - (-1/3)| = {abs(g_abc - r_ref):.6e}"
        )
        print(
            f"[TH-5 1a2 (O2)] by projection: Gamma = {p_abc:.6e}, |Gamma - (-1/3)| = "
            f"{abs(p_abc - r_ref):.6e} (printed; predicted 2.6e-4 x/ 3)"
        )
        print(
            f"[TH-5 1a2 (O2)] lambda0/40 deg 1: |Gamma - (-1/3)| = {abs(g_abc40 - r_ref):.6e}, "
            f"projection {abs(p_abc40 - r_ref):.6e} (printed)"
        )
        print(
            f"[TH-5 1a2 (O2)] lambda0/40 deg 2: |Gamma - (-1/3)| = {abs(g_abc40_d2 - r_ref):.6e}, "
            f"projection {abs(p_abc40_d2 - r_ref):.6e} (printed)"
        )
        g0 = g_pec * np.exp(-2j * BETA_O * L_O)
        print(
            f"[TH-5 1a2 (O) control PEC @ lambda0/80] Gamma = {g_pec:.6f}, "
            f"|1+Gamma(0)| = {abs(1 + g0):.4f}, "
            f"factor |Gamma_PEC - R|/|Gamma_ABC - R| = {abs(g_pec - r_ref) / abs(g_abc - r_ref):.1f}"
        )
    assert abs(g_exact) <= O_BAND, g_exact
    assert abs(g_abc - r_ref) <= O_BAND, (g_abc, r_ref)
