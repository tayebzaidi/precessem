"""`IMP-6` step 2 (§9 item 408; attempt 2 = §9 item 426): the rim resolution
ladder at aspect 5.

Step 1's aspect-5 in-plane fixture (PEC oblate spheroid a = 6 mm, c = 1.2 mm,
ρ_r = c²/a = 0.24 mm, in `IMP-1`'s lossy host, degree 2) solved on three meshes
that differ only in ``resolution_tip`` at ρ_r/3, ρ_r/4 and ρ_r/5 (the rung-2
divisor from ``FEM_EM_IMP6_RUNG2``, the brief's pre-stated 4.5 fallback).
``resolution_surface`` = a/12 and ``resolution_far`` = a/4 are held fixed, so
this is a rim-local ladder, **not** a convergence study of the fixture.

Attempt 1 (item 408) parked: on the default oblate path the finer rungs put
31 / 93 air mask-cell midpoints inside the true spheroid (the rim truncated).
Attempt 2 (item 426, the 2026-10-13 03:00 review ruling) meshes with a
licensed variant of the private oblate builder, chosen by
``FEM_EM_IMP6_RIM`` in {``edge`` (body by revolution, the equator a model
edge), ``tube`` (rim Balls graded on rim distance), ``balls-nonetgen``
(step-1 construction, Netgen pass dropped on rungs ρ_r/4, ρ_r/5)}.
W0 asserts per rung the cell guard, a non-empty mask and **zero** mask-cell
midpoints inside the true spheroid; ``test_census_control`` asserts the
default path still reads > 0 on rungs ρ_r/4 and ρ_r/5 (attempt 1's 31 / 93).

Gated per rung at imported bands (as step 1): (i) |β − 1| ≤ `TH15_BETA_BAND`,
(ii) |Φ_rim,h/Φ_rim,cf,h − 1| ≤ `MAT4_BAND`, (iii) cavity dofs zero.
Printed, not asserted: the signed miss per rung, measured h_rim/ρ_r, and the
inputs for the pairwise observed order (no imported bar backs a rate here).

Step 2b (§9 item 440, chain I26, the 2026-10-13 18:00 daily's provisional
ruling): one rung (ρ_r/4, ``edge``) separating the two candidates for item
426's +1.2–1.6 % floor. (A) the fixed surface grain: ``FEM_EM_IMP6_SURFACE_DIV``
(default 12, bitwise the landed a/12) sets ``resolution_surface`` = a/div and is
in the mesh cache key; ``test_surface_census`` (W0) gates the a/18 mesh before
``test_rung1`` (W1s) solves on it. (B) the comparand's degree-8 quadrature:
``test_quadrature_rung1`` (W1q, no solve) reads Φ_cf,h at degrees 8 / 12 / 14 on
item 426's rung-1 mesh. Asserted: the instrument identities (Q8), (M′_q), (Qh)
and the landed bands; the attribution itself is printed, never asserted.

Step 2c (§9 item 441, chain I27, accepted by the 2026-10-14 03:00 daily ruling):
one rung (ρ_r/4, ``edge``) separating the two remaining candidates for the
floor. (C) the far grain: ``FEM_EM_IMP6_FAR_DIV`` (default 4, bitwise the landed
a/4) sets ``resolution_far`` = a/div and is in the mesh cache key;
``test_surface_census`` (W0f) gates the a/5 mesh ((Kf), cells ≥ 1.05 × 53 959)
before ``test_rung1`` (W1f) solves on it. (D) the quasi-static comparand's
finite-k error: ``FEM_EM_IMP6_K_SCALE`` = s (default 1) scales f and σ
together, so k² → s²k² at a fixed loss tangent; it is **not** in the mesh cache
key (the mesh is k-independent). W1k / W1k2 assert (K) |Im E|/|Re E| ≤
2·s²·8.523e-3 and (Wk) the comparand bitwise on item 440's mesh. One knob per
run (``ValueError`` at import). ``test_knob_defaults`` (W0r, no solve) asserts
the default path bitwise against item 440's records. The (C) / (D) attribution
is printed, never asserted.

Step 2d (§9 item 450, chain I28, the 2026-10-14 10:30 daily's provisional
step): I22's whole rim ladder (ρ_r/3, ρ_r/4, ρ_r/5, ``edge``) re-solved at
``FEM_EM_IMP6_K_SCALE`` = 0.1, so the k → 0 residual is read on every rung.
Under ``K_SCALE != 1`` only: (Wk) stays item 440's 1e-12 identity on rung ρ_r/4;
rungs ρ_r/3 and ρ_r/5 assert (Wk′) |Φ_cf,h(8)/item 426's per-rung record − 1| ≤
5e-8 (the 6-decimal record's half-unit × ≈ 2); rung ρ_r/4 at s = 0.1 asserts (R4)
|Φ_h/22.852453 − 1| ≤ 5e-8 against item 441's W1k (an instrument identity).
``FEM_EM_IMP6_RUNG2`` must stay unset under the knob (``ValueError`` at import).
The pre-registered (Disc) / (Fix) / neither reading rule is in the brief
(docs/planning/chunks/IMP-6.md, "§9 item 450") and is computed from the three
windows' logs, never coded or asserted here. The default path is unchanged.
"""

from __future__ import annotations

import os
import time

import basix
import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    CAVITY_TAG,
    DEGREE,
    E0,
    EPSILON_R_HOST,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    SPHERE_RADIUS,
    TH15_BETA_BAND,
    _assemble_real,
    _cavity_dofs,
)
from tests.validation.test_imp6_oblate_fem import (  # private names only (never test_*)
    BOX_HALF_EXTENTS,
    CELL_GUARD,
    QDEG,
    SEMI_MAJOR,
    _beta,
    _body_volume_deficit,
    _dirichlet,
    _field,
    _imp1_dipole_numpy,
    _ncells,
    _numpy_quadrature_mean,
    _owned_masked_cells,
    _probe_points,
    _rho_r,
    _rim_mask,
    _rss_gib,
    _semi_minor,
)
from tests.validation.test_lossy_sphere_sar import FREQUENCY_HZ, SIGMA_HIGH

ASPECT = 5.0
RIM = os.environ.get("FEM_EM_IMP6_RIM", "edge")
if RIM not in ("edge", "tube", "balls-nonetgen"):
    raise ValueError(f"FEM_EM_IMP6_RIM must be edge, tube or balls-nonetgen, got {RIM!r}")
DIVISORS = (3.0, 4.0, float(os.environ.get("FEM_EM_IMP6_RUNG2", "5")))
# Attempt 1's default-path inside counts (20261013T033815Z_IMP-6-step2-W0diag.log:279, :423).
ATTEMPT1_INSIDE = {4.0: 31, 5.0: 93}
# Step 1's rung-0 record (20261012T184537Z_IMP-6-step1-W1b.log:189-191), printed beside.
STEP1_RUNG0 = dict(cells=32_676, ratio=1.014181, beta=1.006343)
PREDICTED_MISS = {
    3.0: "[0.3 %, 2.0 %] beside step 1's 1.418 %",
    4.0: "~0.8 %",
    5.0: "~0.5 %, range [0.3 %, 1.4 %]",
    4.5: "range [0.3 %, 1.4 %]",
}
# Attempt 1's measured cells (20261013T033117Z_IMP-6-step2-W0.log); item 426: +-15 %.
PREDICTED_CELLS = {3.0: 32_676, 4.0: 51_898, 5.0: 87_426, 4.5: 72_000}

# Item 440 (step 2b): the surface-grain knob (A). float("12") == 12.0, so the
# default is bitwise the landed resolution_surface = SEMI_MAJOR / 12.0.
SURFACE_DIV = float(os.environ.get("FEM_EM_IMP6_SURFACE_DIV", "12"))
if not SURFACE_DIV > 0:
    raise ValueError(f"FEM_EM_IMP6_SURFACE_DIV must be > 0, got {SURFACE_DIV!r}")
# Item 441 (step 2c): the far-grain knob (C) and the k-scale knob (D).
# float("4") == 4.0 and x * 1.0 == x, so both defaults are bitwise the landed path.
FAR_DIV = float(os.environ.get("FEM_EM_IMP6_FAR_DIV", "4"))
if not FAR_DIV > 0:
    raise ValueError(f"FEM_EM_IMP6_FAR_DIV must be > 0, got {FAR_DIV!r}")
K_SCALE = float(os.environ.get("FEM_EM_IMP6_K_SCALE", "1"))
if not 0.0 < K_SCALE <= 1.0:
    raise ValueError(f"FEM_EM_IMP6_K_SCALE must be in (0, 1], got {K_SCALE!r}")
_ACTIVE_KNOBS = [
    n
    for n, on in (
        ("SURFACE_DIV", SURFACE_DIV != 12.0),
        ("FAR_DIV", FAR_DIV != 4.0),
        ("K_SCALE", K_SCALE != 1.0),
    )
    if on
]
if len(_ACTIVE_KNOBS) > 1:
    raise ValueError(f"one knob per run (item 441); active: {_ACTIVE_KNOBS}")
# Item 440's default-path records on rung rho_r/4 (W0r-c):
# 20261014T022824Z_IMP-6-step2b-W0r.log:292; 20261014T022856Z_IMP-6-step2b-W1q.log:177-178.
ITEM440_R1 = dict(cells=53_959, n_mask=11_468, phi_cf_h8=22.704739974621)
# Item 426's W1b |Im E|/|Re E| (20261013T200607Z_IMP-6-step2-a2-W1b.log:178), the (K) basis.
ITEM426_IMAG_RATIO = 8.523e-3
# Step 1's k (20261012T184107Z_IMP-6-step1-W0.log:392), printed beside.
STEP1_K = dict(abs_k=17.900, ka2=0.011535)
# Item 441's pre-registered predictions (rule (e); printed, never asserted).
ITEM441_PRED_MISS_D = "m_D in [0.4 %, 1.6 %]; (D) fires iff m_D <= 0.986 %"
ITEM441_PRED_MISS_C = "m_C in [1.1 %, 1.6 %]; (C) fires iff m_C <= 0.986 %"
ITEM441_PRED_CELLS_C = "~66 k (range 63-68 k; bound ~71-72 k)"
ITEM441_THRESHOLD = 0.00986
ITEM440_M_A = 0.013712
# Item 440 (B): comparand quadrature degrees read; the landed gate keeps QDEG (8).
QDEGS_2B = (8, 12, 14)
# Item 426's rung-1 record (20261013T200607Z_IMP-6-step2-a2-W1b.log:178-179), printed beside.
ITEM426_RUNG1 = dict(phi_h=23.042163, phi_cf_h=22.704740, miss=0.014861, beta=0.992222)
# Item 440's pre-registered predictions (rule (e); printed, never asserted).
ITEM440_PRED_MISS_A = "m_A in [0.5 %, 1.6 %]; (A) fires iff m_A <= 0.986 %"
ITEM440_PRED_CELLS = "~60-65 k (bound ~73-75 k)"
ITEM440_PRED_MB = "m_B <~ 0.05 %; (B) fires iff m_B >= 0.6 %"

# Item 450 (step 2d, I28): item 426's per-rung s = 1 records (printed 6 decimals).
# rho_r/3: 20261013T200520Z_IMP-6-step2-a2-W1a.log:169-171;
# rho_r/4: 20261013T200607Z_IMP-6-step2-a2-W1b.log:177-179;
# rho_r/5: 20261013T200729Z_IMP-6-step2-a2-W1c.log:190-192.
ITEM426_RUNGS = {
    3.0: dict(
        phi_h=22.971085,
        phi_cf_h=22.693485,
        miss=0.012233,
        beta=1.001961,
        imag_ratio=8.521e-3,
        h_max=0.6909,
        h_median=0.5115,
        log="20261013T200520Z_IMP-6-step2-a2-W1a.log:170-171",
    ),
    4.0: dict(
        phi_h=23.042163,
        phi_cf_h=22.704740,
        miss=0.014861,
        beta=0.992222,
        imag_ratio=8.523e-3,
        h_max=0.5103,
        h_median=0.3873,
        log="20261013T200607Z_IMP-6-step2-a2-W1b.log:178-179",
    ),
    5.0: dict(
        phi_h=23.079024,
        phi_cf_h=22.712699,
        miss=0.016129,
        beta=0.994011,
        imag_ratio=8.527e-3,
        h_max=0.3944,
        h_median=0.3089,
        log="20261013T200729Z_IMP-6-step2-a2-W1c.log:191-192",
    ),
}
# (Wk') / (R4) bar: the 6-decimal print's half-unit 5e-7/22.69 = 2.2e-8, times ~2.
ITEM450_REPRO_BAR = 5e-8
# (R4): item 441's W1k Phi_h on rung rho_r/4 at s = 0.1
# (20261014T094547Z_IMP-6-step2c-W1k.log:181).
ITEM441_W1K_PHI_H = 22.852453
# Item 450's pre-registered per-rung predictions at s = 0.1 (rule (e); printed only).
ITEM450_PRED_MISS = {
    3.0: "r3 in [0.34 %, 0.44 %]; delta_D within +-0.05 pp of +0.844 pp",
    4.0: "r4 = +0.6506 % (reproducing W1k, (R4)); delta_D within +-0.05 pp of +0.844 pp",
    5.0: "r5 in [0.73 %, 0.83 %]; delta_D within +-0.05 pp of +0.844 pp",
}
ITEM450_PRED_BETA = {3.0: 0.998, 4.0: 0.988, 5.0: 0.990}
if K_SCALE != 1.0 and DIVISORS[2] != 5.0:
    raise ValueError(
        f"item 450: FEM_EM_IMP6_RUNG2 must stay unset under FEM_EM_IMP6_K_SCALE "
        f"(records exist for rho_r/5 only), got rung-2 divisor {DIVISORS[2]!r}"
    )

_MESH_CACHE: dict = {}


def _variant_kwargs(variant, divisor):
    if variant == "default":
        return {}
    if variant == "balls-nonetgen":
        return dict(oblate_rim="balls", oblate_netgen=(divisor == 3.0))
    return dict(oblate_rim=variant)


def _mesh_rung(divisor, variant=None):
    variant = RIM if variant is None else variant
    key = (variant, divisor, SURFACE_DIV, FAR_DIV)
    if key not in _MESH_CACHE:
        t0 = time.perf_counter()
        out = MeshGenerator.spheroid_in_box_domain(
            semi_major=SEMI_MAJOR,
            semi_minor=_semi_minor(ASPECT),
            box_half_extents=BOX_HALF_EXTENTS,
            resolution_tip=_rho_r(ASPECT) / divisor,
            resolution_surface=SEMI_MAJOR / SURFACE_DIV,
            resolution_far=SEMI_MAJOR / FAR_DIV,
            comm=MPI.COMM_WORLD,
            oblate=True,
            **_variant_kwargs(variant, divisor),
        )
        _MESH_CACHE[key] = (out, time.perf_counter() - t0)
    return _MESH_CACHE[key][0]


def _mesh_seconds(divisor, variant=None):
    return _MESH_CACHE[(RIM if variant is None else variant, divisor, SURFACE_DIV, FAR_DIV)][1]


def _k_numbers():
    """Scaled f, sigma and k = omega sqrt(mu0 eps0 (eps_r - j sigma/(omega eps0))); pure arithmetic."""
    eps0 = 8.8541878128e-12
    mu0 = 1.25663706212e-6
    f = FREQUENCY_HZ * K_SCALE
    sigma = SIGMA_HIGH * K_SCALE
    omega = 2.0 * np.pi * f
    k2 = omega**2 * mu0 * eps0 * (EPSILON_R_HOST - 1j * sigma / (omega * eps0))
    k = complex(np.sqrt(k2))
    return dict(f=f, sigma=sigma, k=k, k2=k2, abs_k=abs(k), ka2=(abs(k) * SEMI_MAJOR) ** 2)


def _numpy_quadrature_mean_q(msh, cells, field_fn, qdeg):
    """Step 1's ``_numpy_quadrature_mean`` with the basix degree as an argument."""
    dmap = msh.geometry.dofmap
    if dmap.shape[1] != 4:
        raise RuntimeError(f"expected affine tets, geometry dofmap {dmap.shape}")
    pts_ref, wts = basix.make_quadrature(basix.CellType.tetrahedron, qdeg)
    xg = msh.geometry.x
    v = xg[dmap[cells]]  # (n, 4, 3)
    jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
    det = np.abs(np.linalg.det(jac)) if len(cells) else np.zeros(0)
    phys = v[:, 0][:, None, :] + np.einsum("cij,qj->cqi", jac, pts_ref)
    if len(cells):
        e = field_fn(phys.reshape(-1, 3)).reshape(len(cells), len(wts), 3)
        e2 = np.sum(np.abs(e) ** 2, axis=2)
        num = float(np.sum(det[:, None] * wts[None, :] * e2))
        vol = float(np.sum(det) * np.sum(wts))
    else:
        num = vol = 0.0
    num = msh.comm.allreduce(num, op=MPI.SUM)
    vol = msh.comm.allreduce(vol, op=MPI.SUM)
    return num / (E0**2 * vol)


def _ufl_dipole_mean_q(msh, chi, qdeg):
    """`IMP-1`'s ``_phi_closed_form_on_mask`` at β = 1 with ``quadrature_degree`` qdeg."""
    comm = msh.comm
    x = ufl.SpatialCoordinate(msh)
    r2 = ufl.dot(x, x)
    r = ufl.sqrt(r2)
    z_hat = ufl.as_vector((0.0, 0.0, 1.0))
    e_cf = E0 * z_hat + E0 * SPHERE_RADIUS**3 * (3.0 * x[2] * x / r**5 - z_hat / r**3)
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": qdeg})
    num = _assemble_real(chi * ufl.dot(e_cf, e_cf) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    return num / (E0**2 * vol)


def _quadrature_identities(msh, chi, cells, label):
    """(M′_q) for q in {12, 14}: numpy vs UFL quadrature of IMP-1's dipole, asserted ≤ 1e-12."""
    comm = msh.comm
    for q in (12, 14):
        mq = _numpy_quadrature_mean_q(msh, cells, _imp1_dipole_numpy, q)
        mu = _ufl_dipole_mean_q(msh, chi, q)
        rel = abs(mq / mu - 1.0)
        if comm.rank == 0:
            print(
                f"\n[IMP-6 step 2b M'_q] {label}: degree {q}: numpy {mq:.15f} vs UFL "
                f"{mu:.15f}, rel {rel:.3e} (asserted <= 1e-12; step 1's degree-8 5.5e-14)",
                flush=True,
            )
        assert rel <= 1e-12, (label, q, mq, mu)


def _inside_count(msh, cells, comm):
    """Owned mask-cell midpoints inside the true spheroid (reduced) and the depth
    a - min rho of the deepest (attempt 1's W0diag computation)."""
    c = _semi_minor(ASPECT)
    mids = (
        msh.geometry.x[msh.geometry.dofmap[cells]].mean(axis=1) if cells.size else np.zeros((0, 3))
    )
    s = (mids[:, 0] ** 2 + mids[:, 1] ** 2) / SEMI_MAJOR**2 + mids[:, 2] ** 2 / c**2
    n_in = int(comm.allreduce(int(np.sum(s < 1.0)), op=MPI.SUM))
    rho_in = np.hypot(mids[:, 0], mids[:, 1])[s < 1.0]
    depth = comm.allreduce(float(SEMI_MAJOR - np.min(rho_in)) if rho_in.size else 0.0, op=MPI.MAX)
    return n_in, depth


def _near_rim_edges(msh, cells):
    """(longest edge, midpoint) of owned mask cells with midpoint within ρ_r of a rim point."""
    rho_r = _rho_r(ASPECT)
    if not cells.size:
        return np.zeros(0), np.zeros((0, 3))
    v = msh.geometry.x[msh.geometry.dofmap[cells]]
    edges = np.max(
        [np.linalg.norm(v[:, i] - v[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)],
        axis=0,
    )
    mid = v.mean(axis=1)
    near = (
        np.minimum(
            np.linalg.norm(mid - [SEMI_MAJOR, 0, 0], axis=1),
            np.linalg.norm(mid + [SEMI_MAJOR, 0, 0], axis=1),
        )
        < rho_r
    )
    return edges[near], mid[near]


def _h_rim_over_rho(msh, cells, comm):
    """Longest edge of mask cells with midpoint within ρ_r of a rim point, / ρ_r (step 1's H)."""
    edges, _ = _near_rim_edges(msh, cells)
    h_local = float(np.max(edges)) if edges.size else 0.0
    return comm.allreduce(h_local, op=MPI.MAX) / _rho_r5()


def _rho_r5():
    return _rho_r(ASPECT)


def _h_rim_stats(msh, cells, comm):
    """Diagnostics beside step 1's max-edge H: median / p95 of the near-rim longest
    edges (/ρ_r), their count, and the midpoint of the max-edge cell (gathered)."""
    edges, mids = _near_rim_edges(msh, cells)
    all_e = np.concatenate(comm.allgather(edges))
    all_m = np.concatenate(comm.allgather(mids.reshape(-1, 3)))
    rho = _rho_r5()
    if not all_e.size:
        return dict(n=0, median=0.0, p95=0.0, where=None)
    return dict(
        n=int(all_e.size),
        median=float(np.median(all_e)) / rho,
        p95=float(np.percentile(all_e, 95)) / rho,
        where=all_m[int(np.argmax(all_e))],
    )


def _solve_on(divisor):
    """Step 1's ``_solve_oblate`` on rung ``divisor``, in-plane drive."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = _mesh_rung(divisor)
    if comm.rank == 0:  # pure arithmetic, nothing collective
        kn = _k_numbers()
        print(
            f"\n[IMP-6 step 2c k] K_SCALE {K_SCALE:g}: f = {kn['f'] / 1e6:.4f} MHz, sigma = "
            f"{kn['sigma']:.5f} S/m, k = {kn['k'].real:.4f} {kn['k'].imag:+.4f}j 1/m, Re k^2 = "
            f"{kn['k2'].real:.4f}, Im k^2 = {kn['k2'].imag:.4f} 1/m^2, |k| = {kn['abs_k']:.4f} 1/m, "
            f"(|k|a)^2 = {kn['ka2']:.4%} (step 1: |k| {STEP1_K['abs_k']:.3f} 1/m, (|k|a)^2 "
            f"{STEP1_K['ka2']:.4%}; 20261012T184107Z_IMP-6-step1-W0.log:392)",
            flush=True,
        )
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ * K_SCALE,
        material=HomogeneousMaterial(sigma=SIGMA_HIGH * K_SCALE, epsilon_r=EPSILON_R_HOST),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=_dirichlet(ASPECT, "inplane"),
        pec_facet_tags=(OUTER_BOUNDARY_TAG, CAVITY_TAG),
    )
    solver = TimeHarmonicSolver(problem, degree=DEGREE)
    fields = solver.solve()

    points = _probe_points(ASPECT)
    re_vals, ok_re = evaluate_vector_field_parallel(fields.e_real, points, comm)
    im_vals, ok_im = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    if not np.all(ok_re & ok_im):
        raise RuntimeError("probe points were not evaluated")
    e_re = np.real(re_vals[:, :3])
    e_im = np.real(im_vals[:, :3])
    # In-plane basis, uniform part on axis 0 (the drive).
    beta = _beta(e_re, points, ASPECT, "inplane", 0)
    imag_ratio = float(
        np.max(np.linalg.norm(e_im, axis=1))
        / max(float(np.max(np.linalg.norm(e_re, axis=1))), 1e-300)
    )

    chi = _rim_mask(msh, ASPECT)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    phi_h = num / (E0**2 * vol)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    phi_cf_h = _numpy_quadrature_mean(msh, cells, _field(ASPECT, "inplane", "continue"))
    # Item 440: comparand at degrees 12 / 14 (printed) and Phi_h at degree 14 ((Qh)).
    phi_cf_q = {
        q: _numpy_quadrature_mean_q(msh, cells, _field(ASPECT, "inplane", "continue"), q)
        for q in (12, 14)
    }
    dx14 = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 14})
    phi_h14 = _assemble_real(chi * ufl.inner(e, e) * dx14, comm) / (
        E0**2 * _assemble_real(chi * dx14, comm)
    )

    cavity = _cavity_dofs(msh, facet_tags, solver.function_space())
    n_cavity = int(comm.allreduce(int(cavity.size), op=MPI.SUM))
    sol = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(sol[cavity]))) if cavity.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))

    return dict(
        beta=beta,
        imag_ratio=imag_ratio,
        phi_h=phi_h,
        phi_cf_h=phi_cf_h,
        phi_cf_q=phi_cf_q,
        phi_h14=phi_h14,
        n_mask=n_mask,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=_ncells(msh),
        h_rim_over_rho=_h_rim_over_rho(msh, cells, comm),
        h_stats=_h_rim_stats(msh, cells, comm),
        deficit=_body_volume_deficit(msh, ASPECT),
        mesh_seconds=_mesh_seconds(divisor),
        seconds=time.perf_counter() - t0,
        rss=_rss_gib(comm),
    )


# --------------------------------------------------------------------- W0
@complex_only
def test_ladder_census():
    comm = MPI.COMM_WORLD
    for divisor in DIVISORS:
        msh, _, _ = _mesh_rung(divisor)
        ncells = _ncells(msh)
        chi = _rim_mask(msh, ASPECT)
        cells = _owned_masked_cells(msh, chi)
        n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
        h = _h_rim_over_rho(msh, cells, comm)
        hs = _h_rim_stats(msh, cells, comm)
        # (G1): mask-cell midpoints inside the true spheroid, i.e. air cells whose
        # faceting truncates the rim (attempt 1: 0 / 31 / 93 on the default path).
        c = _semi_minor(ASPECT)
        n_in, depth = _inside_count(msh, cells, comm)
        if comm.rank == 0:
            print(
                f"\n[IMP-6 step 2 census] variant {RIM}: rung rho_r/{divisor:g} mask-cell midpoints inside the "
                f"true spheroid: {n_in}; deepest at rho = a - {depth * 1e3:.4f} mm "
                f"(focal ring at a - {(SEMI_MAJOR - np.sqrt(SEMI_MAJOR**2 - c**2)) * 1e3:.4f} mm)",
                flush=True,
            )
        deficit = _body_volume_deficit(msh, ASPECT)
        rss = _rss_gib(comm)  # collective: outside the rank-0 branch
        if comm.rank == 0:
            print(
                f"\n[IMP-6 step 2 census] rung rho_r/{divisor:g} near-rim longest edges: n "
                f"{hs['n']}, median/rho_r {hs['median']:.3f}, p95/rho_r {hs['p95']:.3f}, "
                f"max-edge cell midpoint {hs['where']}",
                flush=True,
            )
            print(
                f"\n[IMP-6 step 2 census] rung rho_r/{divisor:g}: resolution_tip "
                f"{_rho_r(ASPECT) / divisor * 1e3:.4f} mm, cells {ncells} (predicted "
                f"~{PREDICTED_CELLS.get(divisor)} +-15 %; guard {CELL_GUARD}), mesh "
                f"{_mesh_seconds(divisor):.1f} s, rim-mask cells {n_mask}, h_rim/rho_r "
                f"{h:.3f}, body-volume deficit {deficit:.3e}, RSS(sum) {rss:.2f} GiB, "
                f"surface divisor {SURFACE_DIV:g} (resolution_surface "
                f"{SEMI_MAJOR / SURFACE_DIV * 1e3:.4f} mm)",
                flush=True,
            )
        assert ncells <= CELL_GUARD, f"cell guard: rung {divisor}: {ncells} > {CELL_GUARD}"
        assert n_mask > 0
        assert n_in == 0, f"(G1) rung rho_r/{divisor:g}, variant {RIM}: {n_in} inside midpoints"


@complex_only
def test_census_control():
    """(C-fid): the default oblate path still truncates the rim on rungs rho_r/4, rho_r/5."""
    comm = MPI.COMM_WORLD
    for divisor in (4.0, 5.0):
        msh, _, _ = _mesh_rung(divisor, "default")
        cells = _owned_masked_cells(msh, _rim_mask(msh, ASPECT))
        n_in, depth = _inside_count(msh, cells, comm)
        ncells = _ncells(msh)
        if comm.rank == 0:
            print(
                f"\n[IMP-6 step 2 control] default path, rung rho_r/{divisor:g}: cells "
                f"{ncells}, inside midpoints {n_in} (attempt 1: {ATTEMPT1_INSIDE[divisor]}), "
                f"deepest at a - {depth * 1e3:.4f} mm",
                flush=True,
            )
        assert n_in > 0, f"(C-fid) default path rung rho_r/{divisor:g} reads 0 inside"


@complex_only
def test_quadrature_rung1():
    """W1q (item 440, (B)): the comparand at degrees 8 / 12 / 14 on rung rho_r/4, no solve."""
    comm = MPI.COMM_WORLD
    msh, _, _ = _mesh_rung(4.0)
    chi = _rim_mask(msh, ASPECT)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    fn = _field(ASPECT, "inplane", "continue")
    phi8_ref = _numpy_quadrature_mean(msh, cells, fn)
    phi = {q: _numpy_quadrature_mean_q(msh, cells, fn, q) for q in QDEGS_2B}
    q8 = abs(phi[8] / phi8_ref - 1.0)
    ncells = _ncells(msh)  # collective: outside the rank-0 branch (first W1q deadlocked here)
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 2b W1q] variant {RIM}, rung rho_r/4, surface divisor "
            f"{SURFACE_DIV:g} (resolution_surface {SEMI_MAJOR / SURFACE_DIV * 1e3:.4f} mm): "
            f"cells {ncells}, rim-mask cells {n_mask}",
            flush=True,
        )
        for q in QDEGS_2B:
            npts = len(basix.make_quadrature(basix.CellType.tetrahedron, q)[1])
            print(
                f"  Phi_cf,h({q}) = {phi[q]:.12f} ({npts} points/cell), "
                f"Phi_cf,h({q})/Phi_cf,h(8) - 1 = {phi[q] / phi[8] - 1:+.3e}",
                flush=True,
            )
        print(
            f"  item 426's degree-8 record {ITEM426_RUNG1['phi_cf_h']:.6f} "
            f"(20261013T200607Z_IMP-6-step2-a2-W1b.log:179; printed beside, not asserted); "
            f"(Q8) |helper(8)/imported - 1| = {q8:.3e} (asserted <= 1e-13)",
            flush=True,
        )
        m_b = abs(phi[14] / phi[8] - 1.0)
        print(
            f"  [step 2b decision (B)] m_B = {m_b:.4%} (degree 12: "
            f"{abs(phi[12] / phi[8] - 1.0):.4%}; predicted {ITEM440_PRED_MB}); (B) explains "
            f"the floor: {'YES' if m_b >= 0.006 else 'NO'} (printed, never asserted)",
            flush=True,
        )
    assert q8 <= 1e-13, ("(Q8)", phi[8], phi8_ref)
    _quadrature_identities(msh, chi, cells, f"rung rho_r/4, surface a/{SURFACE_DIV:g}")


@complex_only
def test_surface_census():
    """W0 (item 440, (A)): rung rho_r/4 only at the module's SURFACE_DIV; no solve."""
    comm = MPI.COMM_WORLD
    msh, _, _ = _mesh_rung(4.0)
    ncells = _ncells(msh)
    chi = _rim_mask(msh, ASPECT)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    h = _h_rim_over_rho(msh, cells, comm)
    hs = _h_rim_stats(msh, cells, comm)
    n_in, depth = _inside_count(msh, cells, comm)
    deficit = _body_volume_deficit(msh, ASPECT)
    phi_cf = {
        q: _numpy_quadrature_mean_q(msh, cells, _field(ASPECT, "inplane", "continue"), q)
        for q in QDEGS_2B
    }
    rss = _rss_gib(comm)  # collective: outside the rank-0 branch
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 2c W0] far divisor {FAR_DIV:g} (resolution_far "
            f"{SEMI_MAJOR / FAR_DIV * 1e3:.4f} mm), surface divisor {SURFACE_DIV:g}"
            + (
                f"; item 441 predicted cells {ITEM441_PRED_CELLS_C}; (Kf) asserted cells >= "
                f"{int(np.ceil(1.05 * ITEM440_R1['cells']))}"
                if FAR_DIV != 4.0
                else ""
            ),
            flush=True,
        )
        print(
            f"\n[IMP-6 step 2b W0] variant {RIM}, rung rho_r/4, surface divisor {SURFACE_DIV:g} "
            f"(resolution_surface {SEMI_MAJOR / SURFACE_DIV * 1e3:.4f} mm): cells {ncells} "
            f"(predicted {ITEM440_PRED_CELLS}; guard {CELL_GUARD}), mesh "
            f"{_mesh_seconds(4.0):.1f} s, rim-mask cells {n_mask}, inside midpoints {n_in} "
            f"(deepest a - {depth * 1e3:.4f} mm), h_rim/rho_r max {h:.3f} / median "
            f"{hs['median']:.3f} / p95 {hs['p95']:.3f}, body-volume deficit {deficit:.3e}, "
            f"RSS(sum) {rss:.2f} GiB (min SICN: the builder's line above)",
            flush=True,
        )
        print(
            f"  comparand on this mesh: Phi_cf,h(8/12/14) = {phi_cf[8]:.12f} / "
            f"{phi_cf[12]:.12f} / {phi_cf[14]:.12f}; m_B here = "
            f"{abs(phi_cf[14] / phi_cf[8] - 1):.4%} (printed beside, outside the rule)",
            flush=True,
        )
    assert ncells <= CELL_GUARD, f"cell guard: {ncells} > {CELL_GUARD}"
    assert n_mask > 0
    assert n_in == 0, f"(G1) a/{SURFACE_DIV:g} rung rho_r/4: {n_in} inside midpoints"
    if FAR_DIV != 4.0:
        assert ncells >= 1.05 * ITEM440_R1["cells"], ("(Kf)", ncells)
    _quadrature_identities(msh, chi, cells, f"rung rho_r/4, surface a/{SURFACE_DIV:g}")


def _predicted_miss(divisor):
    if K_SCALE != 1.0:
        if K_SCALE == 0.1:
            return ITEM450_PRED_MISS.get(divisor, ITEM441_PRED_MISS_D)
        return ITEM441_PRED_MISS_D
    if FAR_DIV != 4.0:
        return ITEM441_PRED_MISS_C
    if SURFACE_DIV != 12.0:
        return ITEM440_PRED_MISS_A
    return PREDICTED_MISS.get(divisor)


@complex_only
def test_knob_defaults():
    """W0r (item 441): the default path re-gated bitwise against item 440's records; no solve."""
    knobs = ("FEM_EM_IMP6_SURFACE_DIV", "FEM_EM_IMP6_FAR_DIV", "FEM_EM_IMP6_K_SCALE")
    set_knobs = [k for k in knobs if k in os.environ]
    if set_knobs:
        pytest.fail(f"(W0r-a) test_knob_defaults cannot run under a knob: {set_knobs}")
    assert SEMI_MAJOR / FAR_DIV == SEMI_MAJOR / 4.0, "(W0r-b) far"
    assert SEMI_MAJOR / SURFACE_DIV == SEMI_MAJOR / 12.0, "(W0r-b) surface"
    assert FREQUENCY_HZ * K_SCALE == FREQUENCY_HZ, "(W0r-b) f"
    assert SIGMA_HIGH * K_SCALE == SIGMA_HIGH, "(W0r-b) sigma"
    comm = MPI.COMM_WORLD
    msh, _, _ = _mesh_rung(4.0)
    ncells = _ncells(msh)
    chi = _rim_mask(msh, ASPECT)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    n_in, _ = _inside_count(msh, cells, comm)
    phi8 = _numpy_quadrature_mean(msh, cells, _field(ASPECT, "inplane", "continue"))
    rel = abs(phi8 / ITEM440_R1["phi_cf_h8"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 2c W0r] default path, rung rho_r/4: cells {ncells} (record "
            f"{ITEM440_R1['cells']}), rim-mask cells {n_mask} (record {ITEM440_R1['n_mask']}), "
            f"inside midpoints {n_in} (record 0), Phi_cf,h(8) = {phi8:.12f} (record "
            f"{ITEM440_R1['phi_cf_h8']:.12f}, rel {rel:.3e}, asserted <= 1e-12)",
            flush=True,
        )
    assert ncells == ITEM440_R1["cells"], ("(W0r-c) cells", ncells)
    assert n_mask == ITEM440_R1["n_mask"], ("(W0r-c) mask", n_mask)
    assert n_in == 0, ("(W0r-c) inside", n_in)
    assert rel <= 1e-12, ("(W0r-c) Phi_cf,h(8)", phi8)


def _report_and_gate(divisor, d):
    comm = MPI.COMM_WORLD
    ratio = d["phi_h"] / d["phi_cf_h"]
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 2] variant {RIM}: rung rho_r/{divisor:g} (aspect 5, in-plane): cells "
            f"{d['ncells']}, mesh {d['mesh_seconds']:.1f} s, solve+post {d['seconds']:.1f} s, "
            f"RSS(sum) {d['rss']:.2f} GiB, rim-mask cells {d['n_mask']}, h_rim/rho_r "
            f"{d['h_rim_over_rho']:.4f} (median {d['h_stats']['median']:.4f}, p95 "
            f"{d['h_stats']['p95']:.4f}), body-volume deficit {d['deficit']:.3e}"
        )
        print(
            f"  beta = {d['beta']:.6f}  |beta - 1| = {abs(d['beta'] - 1):.4%} "
            f"(band {TH15_BETA_BAND:.4%}; predicted within ~0.7 pp of 0.634 %); "
            f"|Im E|/|Re E| = {d['imag_ratio']:.3e}"
        )
        print(
            f"  Phi_rim,h = {d['phi_h']:.6f}, Phi_rim,cf,h = {d['phi_cf_h']:.6f}, "
            f"ratio {ratio:.6f}, signed miss {ratio - 1:+.4%} (band {MAT4_BAND:.0%}; "
            f"predicted "
            f"{_predicted_miss(divisor)})"
        )
        if divisor == 3.0:
            print(
                f"  step-1 rung-0 record (printed beside, not asserted): cells "
                f"{STEP1_RUNG0['cells']}, ratio {STEP1_RUNG0['ratio']:.6f}, beta "
                f"{STEP1_RUNG0['beta']:.6f}"
            )
        qh = abs(d["phi_h14"] / d["phi_h"] - 1.0)
        print(
            f"  [step 2b] surface divisor {SURFACE_DIV:g}, resolution_surface "
            f"{SEMI_MAJOR / SURFACE_DIV * 1e3:.4f} mm; Phi_h(14) = {d['phi_h14']:.12f}, "
            f"|Phi_h(14)/Phi_h(8) - 1| = {qh:.3e} ((Qh), asserted <= 1e-12)"
        )
        for q in (12, 14):
            cfq = d["phi_cf_q"][q]
            print(
                f"  [step 2b] Phi_cf,h({q}) = {cfq:.6f}, Phi_cf,h({q})/Phi_cf,h(8) - 1 = "
                f"{cfq / d['phi_cf_h'] - 1:+.3e}, signed miss vs degree-{q} comparand "
                f"{d['phi_h'] / cfq - 1:+.4%}"
            )
        if SURFACE_DIV != 12.0:
            m_a = ratio - 1.0
            print(
                f"  [step 2b] item 426 rung-1 record (a/12, printed beside): Phi_h "
                f"{ITEM426_RUNG1['phi_h']:.6f}, Phi_cf,h {ITEM426_RUNG1['phi_cf_h']:.6f}, miss "
                f"{ITEM426_RUNG1['miss']:+.4%}, beta {ITEM426_RUNG1['beta']:.6f}"
            )
            print(
                f"  [step 2b decision (A)] m_A = {m_a:+.4%} (predicted {ITEM440_PRED_MISS_A}); "
                f"(A) explains the floor: {'YES' if m_a <= 0.00986 else 'NO'} "
                f"(printed, never asserted)"
            )
        if (K_SCALE != 1.0 and divisor == 4.0) or FAR_DIV != 4.0:
            print(
                f"  [step 2c] item 426 rung-1 record (printed beside): Phi_h "
                f"{ITEM426_RUNG1['phi_h']:.6f}, Phi_cf,h {ITEM426_RUNG1['phi_cf_h']:.6f}, miss "
                f"{ITEM426_RUNG1['miss']:+.4%}, beta {ITEM426_RUNG1['beta']:.6f}, |Im E|/|Re E| "
                f"{ITEM426_IMAG_RATIO:.3e} (20261013T200607Z_IMP-6-step2-a2-W1b.log:178-179)"
            )
        if K_SCALE != 1.0 and divisor != 4.0:
            rec = ITEM426_RUNGS[divisor]
            print(
                f"  [step 2d] item 426 rung rho_r/{divisor:g} record (printed beside): Phi_h "
                f"{rec['phi_h']:.6f}, Phi_cf,h {rec['phi_cf_h']:.6f}, miss {rec['miss']:+.4%}, "
                f"beta {rec['beta']:.6f}, |Im E|/|Re E| {rec['imag_ratio']:.3e} ({rec['log']})"
            )
            s2 = K_SCALE**2
            print(
                f"  [step 2c (K)] |Im E|/|Re E| = {d['imag_ratio']:.3e} (predicted "
                f"{s2 * rec['imag_ratio']:.2e}; asserted <= {2 * s2 * ITEM426_IMAG_RATIO:.3e}); "
                f"(Wk') |Phi_cf,h(8)/{rec['phi_cf_h']:.6f} - 1| = "
                f"{abs(d['phi_cf_h'] / rec['phi_cf_h'] - 1):.3e} (asserted <= "
                f"{ITEM450_REPRO_BAR:.0e})"
            )
        if K_SCALE != 1.0 and divisor in ITEM426_RUNGS:
            rec = ITEM426_RUNGS[divisor]
            s2 = K_SCALE**2
            r_d = ratio - 1.0
            delta_dd = (rec["miss"] - r_d) / (1.0 - s2)
            if K_SCALE == 0.1:
                rule = (
                    f"I28 prediction {ITEM450_PRED_MISS[divisor]}; beta predicted "
                    f"~{ITEM450_PRED_BETA[divisor]:.3f}, |beta - 1| <= 1.4 %"
                )
            else:
                rule = "outside the I28 rule (s != 0.1)"
            print(
                f"  [step 2d] rung rho_r/{divisor:g}, s = {K_SCALE:g}: s = 1 record miss "
                f"{rec['miss']:+.4%} ({rec['log']}); r_d = {r_d:+.4%}; delta_D,d = "
                f"(m_d - r_d)/(1 - s^2) = {delta_dd * 100:+.4f} pp; beta {d['beta']:.6f}; "
                f"{rule}; h_rim/rho_r max {d['h_rim_over_rho']:.4f} / median "
                f"{d['h_stats']['median']:.4f} (s = 1 record {rec['h_max']:.4f} / "
                f"{rec['h_median']:.4f}) (printed, never asserted)"
            )
            if divisor == 4.0 and K_SCALE == 0.1:
                print(
                    f"  [step 2d (R4)] Phi_h = {d['phi_h']:.12f} vs item 441's W1k "
                    f"{ITEM441_W1K_PHI_H:.6f} (20261014T094547Z_IMP-6-step2c-W1k.log:181): "
                    f"|Phi_h/W1k - 1| = {abs(d['phi_h'] / ITEM441_W1K_PHI_H - 1):.3e} "
                    f"(asserted <= {ITEM450_REPRO_BAR:.0e})"
                )
        if K_SCALE != 1.0 and divisor == 4.0:
            m = ratio - 1.0
            s2 = K_SCALE**2
            delta_d = (ITEM426_RUNG1["miss"] - m) / (1.0 - s2)
            print(
                f"  [step 2c (K)] |Im E|/|Re E| = {d['imag_ratio']:.3e} (predicted "
                f"{s2 * ITEM426_IMAG_RATIO:.2e}; asserted <= {2 * s2 * ITEM426_IMAG_RATIO:.3e}); "
                f"(Wk) |Phi_cf,h(8)/{ITEM440_R1['phi_cf_h8']:.12f} - 1| = "
                f"{abs(d['phi_cf_h'] / ITEM440_R1['phi_cf_h8'] - 1):.3e} (asserted <= 1e-12)"
            )
            if K_SCALE == 0.1:
                verdict = "YES" if m <= ITEM441_THRESHOLD else "NO"
                if m - ITEM426_RUNG1["miss"] > 0.005:
                    verdict += " (moves the wrong way)"
            else:
                verdict = "outside the rule (s != 0.1)"
            print(
                f"  [step 2c decision (D)] s = {K_SCALE:g}: m = {m:+.4%} (predicted "
                f"{ITEM441_PRED_MISS_D}); delta_D = (m0 - m)/(1 - s^2) = {delta_d:+.4%}; "
                f"(D) explains the floor: {verdict} (printed, never asserted)"
            )
        if FAR_DIV != 4.0:
            m = ratio - 1.0
            verdict = "YES" if m <= ITEM441_THRESHOLD else "NO"
            if m - ITEM426_RUNG1["miss"] > 0.005:
                verdict += " (moves the wrong way)"
            print(
                f"  [step 2c decision (C)] resolution_far a/{FAR_DIV:g}: m_C = {m:+.4%} (predicted "
                f"{ITEM441_PRED_MISS_C}; item 440's m_A {ITEM440_M_A:+.4%} beside); (C) explains "
                f"the floor: {verdict} (printed, never asserted)"
            )
        print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}", flush=True)
    assert abs(d["phi_h14"] / d["phi_h"] - 1.0) <= 1e-12, ("(Qh)", d["phi_h"], d["phi_h14"])
    if K_SCALE != 1.0:
        assert d["imag_ratio"] <= 2.0 * K_SCALE**2 * ITEM426_IMAG_RATIO, ("(K)", d["imag_ratio"])
        if divisor == 4.0:
            assert abs(d["phi_cf_h"] / ITEM440_R1["phi_cf_h8"] - 1.0) <= 1e-12, (
                "(Wk)",
                d["phi_cf_h"],
            )
            if K_SCALE == 0.1:
                assert abs(d["phi_h"] / ITEM441_W1K_PHI_H - 1.0) <= ITEM450_REPRO_BAR, (
                    "(R4)",
                    d["phi_h"],
                )
        else:
            rec_cf = ITEM426_RUNGS[divisor]["phi_cf_h"]
            assert abs(d["phi_cf_h"] / rec_cf - 1.0) <= ITEM450_REPRO_BAR, (
                "(Wk')",
                divisor,
                d["phi_cf_h"],
            )
    assert d["n_cavity"] > 0
    assert d["cavity_max"] < 1e-12
    assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (divisor, d["beta"], d["h_rim_over_rho"])
    assert abs(ratio - 1.0) <= MAT4_BAND, (divisor, d["phi_h"], d["phi_cf_h"], d["h_rim_over_rho"])


# --------------------------------------------------------------------- W1a
@complex_only
@pytest.mark.integration
def test_rung0():
    _report_and_gate(3.0, _solve_on(3.0))


# --------------------------------------------------------------------- W1b
@complex_only
@pytest.mark.integration
def test_rung1():
    _report_and_gate(4.0, _solve_on(4.0))


# --------------------------------------------------------------------- W1c
@complex_only
@pytest.mark.integration
def test_rung2():
    _report_and_gate(DIVISORS[2], _solve_on(DIVISORS[2]))
