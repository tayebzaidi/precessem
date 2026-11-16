"""PEC rectangular-cavity eigenmodes (TH-9).

A lossless, source-free cavity bounded by a perfect electric conductor is the
purest available check of the curl-curl and mass assembly: the eigenproblem

    ∫ (∇×E)·(∇×v) dx = k² ∫ E·v dx,   n × E = 0 on ∂Ω

is real symmetric — no materials, no sources, no complex scalar type — and its
eigenvalues have a closed form on a box of edges ``a × b × d``:

    f_mnp = (c/2)·√((m/a)² + (n/b)² + (p/d)²)

with at most one of ``m, n, p`` equal to zero. Any sign or scaling error in
either bilinear form moves every eigenvalue, so this gate fails loudly.

Two implementation points that are easy to get wrong:

* **Gradient (null-space) modes.** ``∇×∇φ = 0``, so every discrete gradient is
  an eigenvector with ``k² = 0``. They are physical for the discretisation, not
  numerical noise, and they must be *counted and discarded* against a stated
  cutoff rather than silently skipped — the count is a diagnostic of the N1curl
  null space. See :func:`smallest_cavity_eigenvalues`.
* **Constrained dofs.** The PEC rows are eliminated by assembling ``A`` with a
  large diagonal (``bc_diagonal``, far above the spectrum of interest) and
  ``B`` with unit diagonal. That keeps ``B`` symmetric positive definite — so
  the problem stays a true GHEP — and parks the spurious eigenvalues at
  ``bc_diagonal``, where the cutoff in :func:`solve_pec_cavity_modes` drops
  them.

The shift-and-invert target is derived from the analytic spectrum. That is a
preconditioning choice, not an input to the answer: the returned eigenvalues
are the Rayleigh quotients of the assembled pencil and are compared against the
closed form independently.
"""

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from slepc4py import SLEPc

import dolfinx
from dolfinx import fem, mesh as dmesh
from dolfinx.fem.petsc import assemble_matrix

from ..utils.constants import C_0, EPSILON_0, ETA_0, MU_0


def analytic_cavity_frequencies(
    edges: Sequence[float],
    n_modes: int,
    max_index: int = 4,
) -> np.ndarray:
    """Closed-form PEC rectangular-cavity resonances, ascending [Hz].

    ``edges`` is ``(a, b, d)``. Modes with two or more zero indices do not
    exist (the field vanishes identically), so they are excluded; degenerate
    modes are returned with their multiplicity, because the solver returns them
    that way too.
    """

    return analytic_box_cavity_modes(edges, (), n_modes=n_modes, max_index=max_index)


def frequency_from_eigenvalue(eigenvalue: float) -> float:
    """``k² -> f`` [Hz]. ``k = ω/c``, so ``f = c√(k²)/(2π)``."""

    return C_0 * np.sqrt(max(eigenvalue, 0.0)) / (2.0 * np.pi)


@dataclass
class CavitySpectrum:
    """Result of one cavity eigen-solve."""

    frequencies_hz: np.ndarray
    eigenvalues: np.ndarray
    null_mode_count: int
    null_cutoff: float
    n_converged: int
    n_dofs: int
    n_cells: int
    cell_size: float
    degree: int
    diagnostics: dict = field(default_factory=dict)
    #: Populated only when ``solve_pec_cavity_modes(return_modes=True)``: the
    #: N1curl eigenfunctions of ``frequencies_hz``, in the same order, and the
    #: mesh they live on. Left ``None`` otherwise so the gate path is unchanged.
    mode_functions: list = None
    mesh: object = None


def _cavity_forms(
    V,
    bc_diagonal: float,
    surface_impedance_ohm: complex = None,
    omega_rad_s: float = None,
    pec_facets: np.ndarray = None,
):
    """Assemble the stiffness/mass pencil.

    Default (``surface_impedance_ohm=None``): PEC tangential-E constraints by
    the large-diagonal elimination — this path is unchanged by `TH-14`.

    ``pec_facets`` (`TH-16` step 1): the local facet indices to pin. ``None``
    ⇒ every exterior facet (today's path, byte-identical). A subset leaves the
    remaining exterior facets *natural*, which in the E form is the perfect
    magnetic conductor n × ∇×E = 0 (Jin §8.7.5) — no ``ds`` term. The pin is
    built from the PEC facet list, so an edge shared by a PEC and a PMC face
    (the rim) stays pinned.

    With a complex ``surface_impedance_ohm`` Z_s (`TH-14` step 1), no Dirichlet
    pin is applied; the Leontovich wall ``E_t = Z_s H × n`` (n outward) enters
    as the third-kind boundary term of Jin §5.8.3 on every exterior facet:

        A += jω₀μ₀/Z_s ∫_Γ (n × u)·(n × v) dS,

    from ``∮ (n × ∇×E)·v̄`` with ``∇×E = −jωμ₀H`` (e^{jωt}). The pencil is then
    complex and non-Hermitian; the constrained-DOF count returned is 0.
    """

    msh = V.mesh
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    mass = fem.form(ufl.inner(u, v) * ufl.dx)

    if surface_impedance_ohm is not None:
        if not np.issubdtype(PETSc.ScalarType, np.complexfloating):
            raise RuntimeError(
                "a surface-impedance cavity needs the complex DolfinX build "
                "(source /usr/local/bin/dolfinx-complex-mode)"
            )
        if omega_rad_s is None:
            raise ValueError("surface_impedance_ohm requires omega_rad_s")
        gamma = 1j * omega_rad_s * MU_0 / complex(surface_impedance_ohm)
        gamma_c = fem.Constant(msh, PETSc.ScalarType(gamma))
        n = ufl.FacetNormal(msh)
        stiffness = fem.form(
            ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx
            + gamma_c * ufl.inner(ufl.cross(n, u), ufl.cross(n, v)) * ufl.ds
        )
        A = assemble_matrix(stiffness)
        A.assemble()
        B = assemble_matrix(mass)
        B.assemble()
        return A, B, 0

    stiffness = fem.form(ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx)

    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    if pec_facets is None:
        boundary_facets = dmesh.exterior_facet_indices(msh.topology)
    else:
        boundary_facets = np.asarray(pec_facets, dtype=np.int32)
    boundary_dofs = fem.locate_dofs_topological(V, tdim - 1, boundary_facets)
    zero = fem.Function(V)
    zero.x.array[:] = 0.0
    bc = fem.dirichletbc(zero, boundary_dofs)

    # dolfinx 0.11 renamed this keyword ``diagonal`` -> ``diag``; the semantics
    # are unchanged (BC rows/columns zeroed, the diagonal set to the value), so
    # the constrained-DOF eigenvalues still land at ``bc_diagonal``/1.0 and the
    # ``spurious_cutoff`` reasoning above holds verbatim (`OPS-24`, 2026-08-25).
    A = assemble_matrix(stiffness, bcs=[bc], diag=bc_diagonal)
    A.assemble()
    B = assemble_matrix(mass, bcs=[bc], diag=1.0)
    B.assemble()
    return A, B, boundary_dofs.size


def _solve_pencil(
    A, B, target: float, nev: int, comm, return_vectors: bool = False
) -> Tuple[np.ndarray, int]:
    """Shift-and-invert GHEP solve; returns (ascending eigenvalues, n_converged).

    With ``return_vectors`` the eigenvectors are returned as a third element, in
    the same ascending order as the eigenvalues (SLEPc's own ordering follows the
    shift-and-invert target, not the spectrum). Nothing else changes: no caller
    that leaves the flag alone can observe the difference.
    """

    eps = SLEPc.EPS().create(comm)
    eps.setOperators(A, B)
    eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setTarget(target)
    eps.setDimensions(nev=nev, ncv=max(4 * nev, 20))
    eps.setTolerances(tol=1e-10, max_it=200)

    st = eps.getST()
    st.setType(SLEPc.ST.Type.SINVERT)
    ksp = st.getKSP()
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")

    eps.setFromOptions()
    eps.solve()

    n_converged = eps.getConverged()
    values = np.array(
        [np.real(eps.getEigenvalue(i)) for i in range(n_converged)], dtype=float
    )
    order = np.argsort(values)

    vectors = None
    if return_vectors:
        vectors = []
        for i in order:
            vec = A.createVecRight()
            eps.getEigenvector(int(i), vec)
            vectors.append(vec)

    eps.destroy()
    if return_vectors:
        return values[order], n_converged, vectors
    return values[order], n_converged


def solve_pec_cavity_modes(
    edges: Sequence[float] = (1.0, 0.8, 0.6),
    divisions: Sequence[int] = (8, 7, 5),
    degree: int = 2,
    n_modes: int = 4,
    comm: MPI.Comm = None,
    null_cutoff_fraction: float = 1e-4,
    return_modes: bool = False,
) -> CavitySpectrum:
    """Solve for the lowest ``n_modes`` physical resonances of a PEC box.

    ``null_cutoff_fraction`` is expressed as a fraction of the analytic
    fundamental ``k²``: eigenvalues below it are gradient (null-space) modes and
    are counted, reported, and discarded.

    ``return_modes`` additionally fills ``CavitySpectrum.mode_functions`` and
    ``.mesh`` with the eigenfunctions behind the returned frequencies, so a
    caller can export the very field its assertions were read from. It is purely
    additive: with the flag off, the returned frequencies and every diagnostic
    are bit-identical to before.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    a, b, d = edges
    nx, ny, nz = divisions

    msh = dmesh.create_box(
        comm,
        [np.array([0.0, 0.0, 0.0]), np.array([a, b, d])],
        [nx, ny, nz],
        cell_type=dmesh.CellType.tetrahedron,
    )
    V = fem.functionspace(msh, ("N1curl", degree))

    analytic = analytic_cavity_frequencies(edges, n_modes)
    lam = (2.0 * np.pi * analytic / C_0) ** 2

    # Park the constrained rows well above the band of interest, and target the
    # middle of that band so every requested physical mode is closer to the
    # shift than the gradient cluster at zero is.
    bc_diagonal = 1.0e4 * lam[-1]
    target = 0.5 * (lam[0] + lam[-1])

    A, B, n_bc = _cavity_forms(V, bc_diagonal)
    if return_modes:
        values, n_converged, vectors = _solve_pencil(
            A, B, target, nev=n_modes + 4, comm=comm, return_vectors=True
        )
    else:
        values, n_converged = _solve_pencil(A, B, target, nev=n_modes + 4, comm=comm)

    null_cutoff = null_cutoff_fraction * lam[0]
    spurious_cutoff = 0.5 * bc_diagonal
    null_modes = values[values < null_cutoff]
    keep = (values >= null_cutoff) & (values < spurious_cutoff)
    physical = values[keep]

    if physical.size < n_modes:
        raise RuntimeError(
            f"only {physical.size} physical eigenvalues converged in the band "
            f"[{null_cutoff:.3e}, {spurious_cutoff:.3e}); requested {n_modes}"
        )
    physical = physical[:n_modes]

    mode_functions = None
    if return_modes:
        # Same mask, same order as the frequencies: mode_functions[i] is the
        # eigenfunction of frequencies_hz[i]. The eigenvector's local block is
        # owned-dofs-only, so the ghost values are filled by scatter_forward.
        kept = [vec for vec, take in zip(vectors, keep) if take][:n_modes]
        mode_functions = []
        for index, vec in enumerate(kept):
            f = fem.Function(V, name=f"E_mode_{index + 1}")
            local = f.x.index_map.size_local * f.x.block_size
            f.x.array[:local] = vec.getArray(readonly=True)
            f.x.scatter_forward()
            mode_functions.append(f)
        for vec in vectors:
            vec.destroy()

    n_dofs = V.dofmap.index_map.size_global * V.dofmap.index_map_bs
    n_cells = msh.topology.index_map(msh.topology.dim).size_global
    cell_size = float(max(a / nx, b / ny, d / nz))

    A.destroy()
    B.destroy()

    return CavitySpectrum(
        frequencies_hz=np.array([frequency_from_eigenvalue(v) for v in physical]),
        eigenvalues=physical,
        null_mode_count=int(null_modes.size),
        null_cutoff=float(null_cutoff),
        n_converged=int(n_converged),
        n_dofs=int(n_dofs),
        n_cells=int(n_cells),
        cell_size=cell_size,
        degree=degree,
        diagnostics={
            "analytic_frequencies_hz": analytic,
            "target_eigenvalue": float(target),
            "bc_diagonal": float(bc_diagonal),
            "n_constrained_dofs_local": int(n_bc),
        },
        mode_functions=mode_functions,
        mesh=msh if return_modes else None,
    )


# ---------------------------------------------------------------------------
# `TH-16` step 1: per-face PEC / PMC on the box eigen path.
#
# Additive. ``solve_pec_cavity_modes`` above is untouched; the new solver
# passes an explicit PEC facet list to ``_cavity_forms`` and leaves the faces
# named in ``pmc_faces`` natural (PMC, Jin §8.7.5). The faces at the origin are
# always PEC; a PMC face is one of the three far faces, at most one per axis.
# ---------------------------------------------------------------------------

_PMC_FACE_AXIS = {"x+": 0, "y+": 1, "z+": 2}


def analytic_box_cavity_modes(
    edges: Sequence[float],
    pmc_faces: Sequence[str] = (),
    n_modes: int = 4,
    max_index: int = 6,
) -> np.ndarray:
    """Closed-form resonances of a box with PMC on ``pmc_faces``, ascending [Hz].

    Per axis of length L: ``k = mπ/L`` (PEC–PEC, m ≥ 0) or ``(m + ½)π/L``
    (PEC at 0, PMC at L; m ≥ 0, never zero). The divergence-free amplitude
    space of ``(k_x, k_y, k_z)`` has dimension 2 when all three are non-zero
    and 1 when exactly one is zero (two zeros: the field vanishes), so each
    triple is returned with that multiplicity.
    """

    axes = set()
    for face in pmc_faces:
        if face not in _PMC_FACE_AXIS:
            raise ValueError(f"unknown PMC face {face!r}; expected x+, y+ or z+")
        axis = _PMC_FACE_AXIS[face]
        if axis in axes:
            raise ValueError(f"more than one PMC face on axis {axis}")
        axes.add(axis)

    per_axis = []
    for axis, length in enumerate(edges):
        shift = 0.5 if axis in axes else 0.0
        per_axis.append(
            [(m + shift) * np.pi / length for m in range(max_index + 1)]
        )

    values = []
    for kx in per_axis[0]:
        for ky in per_axis[1]:
            for kz in per_axis[2]:
                zeros = [kx, ky, kz].count(0.0)
                multiplicity = {0: 2, 1: 1}.get(zeros, 0)
                f = C_0 * np.sqrt(kx**2 + ky**2 + kz**2) / (2.0 * np.pi)
                values.extend([f] * multiplicity)
    values = np.sort(np.asarray(values))
    if values.size < n_modes:
        raise ValueError(
            f"max_index={max_index} yields only {values.size} modes; "
            f"{n_modes} were requested"
        )
    return values[:n_modes]


def pmc_face_facets(msh, edges: Sequence[float], pmc_faces: Sequence[str] = ()):
    """Local (rank-local, owned + ghost) facet indices on each named PMC face."""

    tdim = msh.topology.dim
    facets = []
    for face in pmc_faces:
        axis = _PMC_FACE_AXIS[face]
        coord = float(edges[axis])
        facets.append(
            dmesh.locate_entities_boundary(
                msh, tdim - 1, lambda x, a=axis, c=coord: np.isclose(x[a], c)
            )
        )
    if not facets:
        return np.zeros(0, dtype=np.int32)
    return np.unique(np.concatenate(facets)).astype(np.int32)


def solve_box_cavity_modes(
    edges: Sequence[float] = (1.0, 0.8, 0.6),
    divisions: Sequence[int] = (6, 5, 4),
    degree: int = 2,
    n_modes: int = 4,
    pmc_faces: Sequence[str] = (),
    comm: MPI.Comm = None,
    null_cutoff_fraction: float = 1e-4,
) -> CavitySpectrum:
    """Lowest ``n_modes`` resonances of a box with PMC on ``pmc_faces`` (`TH-16`).

    PEC facets = exterior facets minus those on the PMC planes; the PMC faces
    are left natural. Same pencil, target and cutoffs as
    :func:`solve_pec_cavity_modes`, with the shift-invert target taken from
    :func:`analytic_box_cavity_modes`. ``diagnostics["n_pmc_facets_owned"]``
    is the *global* (allreduced) count of owned PMC facets.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    a, b, d = edges
    nx, ny, nz = divisions

    msh = dmesh.create_box(
        comm,
        [np.array([0.0, 0.0, 0.0]), np.array([a, b, d])],
        [nx, ny, nz],
        cell_type=dmesh.CellType.tetrahedron,
    )
    V = fem.functionspace(msh, ("N1curl", degree))

    analytic = analytic_box_cavity_modes(edges, pmc_faces, n_modes)
    lam = (2.0 * np.pi * analytic / C_0) ** 2
    bc_diagonal = 1.0e4 * lam[-1]
    target = 0.5 * (lam[0] + lam[-1])

    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    exterior = dmesh.exterior_facet_indices(msh.topology)
    pmc = pmc_face_facets(msh, edges, pmc_faces)
    pec = np.setdiff1d(exterior, pmc).astype(np.int32)
    n_owned_facets = msh.topology.index_map(tdim - 1).size_local
    n_pmc_owned = comm.allreduce(int(np.sum(pmc < n_owned_facets)), op=MPI.SUM)

    A, B, n_bc = _cavity_forms(V, bc_diagonal, pec_facets=pec)
    values, n_converged = _solve_pencil(A, B, target, nev=n_modes + 4, comm=comm)

    null_cutoff = null_cutoff_fraction * lam[0]
    spurious_cutoff = 0.5 * bc_diagonal
    null_modes = values[values < null_cutoff]
    keep = (values >= null_cutoff) & (values < spurious_cutoff)
    physical = values[keep]
    if physical.size < n_modes:
        raise RuntimeError(
            f"only {physical.size} physical eigenvalues converged in the band "
            f"[{null_cutoff:.3e}, {spurious_cutoff:.3e}); requested {n_modes}"
        )
    physical = physical[:n_modes]

    n_dofs = V.dofmap.index_map.size_global * V.dofmap.index_map_bs
    n_cells = msh.topology.index_map(msh.topology.dim).size_global
    A.destroy()
    B.destroy()

    return CavitySpectrum(
        frequencies_hz=np.array([frequency_from_eigenvalue(v) for v in physical]),
        eigenvalues=physical,
        null_mode_count=int(null_modes.size),
        null_cutoff=float(null_cutoff),
        n_converged=int(n_converged),
        n_dofs=int(n_dofs),
        n_cells=int(n_cells),
        cell_size=float(max(a / nx, b / ny, d / nz)),
        degree=degree,
        diagnostics={
            "analytic_frequencies_hz": analytic,
            "target_eigenvalue": float(target),
            "bc_diagonal": float(bc_diagonal),
            "n_constrained_dofs_local": int(n_bc),
            "pmc_faces": tuple(pmc_faces),
            "n_pmc_facets_owned": int(n_pmc_owned),
            "null_max_abs_over_lam0": (
                float(np.abs(null_modes).max() / lam[0]) if null_modes.size else 0.0
            ),
        },
    )


def _solve_pencil_nonhermitian(
    A, B, target: float, nev: int, comm, return_vectors: bool = False
) -> Tuple[np.ndarray, int]:
    """Shift-and-invert GNHEP solve; complex eigenvalues ordered by |λ − target|.

    `TH-14`: a Hermitian problem type would silently drop Im λ, so the
    lossy-wall pencil (and its PEC control) goes through ``GNHEP``.
    With ``return_vectors`` the eigenvectors follow as a third element, in the
    same |λ − target| order (`EX-60`; additive, as in ``_solve_pencil``).
    """

    eps = SLEPc.EPS().create(comm)
    eps.setOperators(A, B)
    eps.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setTarget(target)
    eps.setDimensions(nev=nev, ncv=max(4 * nev, 20))
    eps.setTolerances(tol=1e-10, max_it=200)
    st = eps.getST()
    st.setType(SLEPc.ST.Type.SINVERT)
    ksp = st.getKSP()
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    eps.setFromOptions()
    eps.solve()
    n_converged = eps.getConverged()
    values = np.array(
        [complex(eps.getEigenvalue(i)) for i in range(n_converged)], dtype=complex
    )
    order = np.argsort(np.abs(values - target))
    vectors = None
    if return_vectors:
        vectors = []
        for i in order:
            vec = A.createVecRight()
            eps.getEigenvector(int(i), vec)
            vectors.append(vec)
    eps.destroy()
    if return_vectors:
        return values[order], n_converged, vectors
    return values[order], n_converged


def surface_resistance_ohm(omega_rad_s: float, sigma_s_per_m: float) -> float:
    """Good-conductor ``R_s = √(ωμ₀/(2σ))`` [Ω]; ``Z_s = (1 + j)R_s``."""

    return float(np.sqrt(omega_rad_s * MU_0 / (2.0 * sigma_s_per_m)))


def pozar_te10l_q_conductor(edges: Sequence[float], ell: int, sigma_s_per_m: float) -> float:
    """Pozar (Microwave Engineering, rectangular-cavity Q_c) for TE₁₀ℓ.

    ``edges = (a, b, d)`` with E along ``b``; R_s taken at the PEC frequency:

        Q_c = (kad)³ b η / (2π² R_s) · 1/(2ℓ²a³b + 2bd³ + ℓ²a³d + ad³)
    """

    a, b, d = edges
    k = np.pi * np.sqrt((1.0 / a) ** 2 + (ell / d) ** 2)
    r_s = surface_resistance_ohm(k * C_0, sigma_s_per_m)
    return float(
        (k * a * d) ** 3 * b * ETA_0 / (2.0 * np.pi**2 * r_s)
        / (2 * ell**2 * a**3 * b + 2 * b * d**3 + ell**2 * a**3 * d + a * d**3)
    )


@dataclass
class ImpedanceWallMode:
    """One complex cavity eigenpair of the lossy-wall (or PEC control) pencil."""

    eigenvalue: complex
    omega_rad_s: complex
    q: float
    omega0_rad_s: float
    n_converged: int
    n_cells: int
    n_dofs: int
    n_constrained_dofs_local: int
    surface_impedance_ohm: complex = None
    mode_function: object = None
    mesh: object = None


def solve_impedance_wall_cavity_mode(
    edges: Sequence[float] = (1.0, 0.8, 0.6),
    divisions: Sequence[int] = (9, 7, 6),
    degree: int = 2,
    sigma_s_per_m: float = None,
    omega_linearisation_rad_s: float = None,
    mode_indices: Tuple[int, int, int] = (1, 0, 1),
    nev: int = 4,
    comm: MPI.Comm = None,
    return_mode: bool = False,
    surface_impedance_fn=None,
) -> ImpedanceWallMode:
    """Complex eigenpair nearest the PEC ``k²`` of ``mode_indices`` (`TH-14`).

    ``sigma_s_per_m=None`` is the negative control: `TH-9`'s own Dirichlet PEC
    pencil, solved non-Hermitian so any spurious ``Im λ`` would be visible.
    Otherwise the Leontovich wall with ``Z_s = (1 + j)R_s(ω_lin)``; ``ω_lin``
    defaults to the PEC ω₀ (the linearisation), and a fixed-point update is a
    second call with ``omega_linearisation_rad_s = Re ω``.
    ``Q = Re ω / (2|Im ω|)`` with ``ω = c√λ``.

    ``surface_impedance_fn`` (`TH-18` step 1): a callable ``ω_lin → Z_s``
    replacing ``(1 + j)R_s(ω_lin)`` on the wall — e.g. a layered stack from
    ``core.layered_impedance``. It is called once, at ``ω_lin``, and takes the
    lossy-wall branch even when ``sigma_s_per_m`` is ``None``. Unset, the path
    is byte-identical to `TH-14`'s.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    a, b, d = edges
    nx, ny, nz = divisions
    msh = dmesh.create_box(
        comm,
        [np.array([0.0, 0.0, 0.0]), np.array([a, b, d])],
        [nx, ny, nz],
        cell_type=dmesh.CellType.tetrahedron,
    )
    V = fem.functionspace(msh, ("N1curl", degree))

    m, n_, p = mode_indices
    k0_sq = float(np.pi**2 * ((m / a) ** 2 + (n_ / b) ** 2 + (p / d) ** 2))
    omega0 = float(np.sqrt(k0_sq) * C_0)

    z_s = None
    if sigma_s_per_m is None and surface_impedance_fn is None:
        A, B, n_bc = _cavity_forms(V, 1.0e4 * k0_sq)
    else:
        omega_lin = omega0 if omega_linearisation_rad_s is None else float(
            omega_linearisation_rad_s
        )
        if surface_impedance_fn is None:
            z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega_lin, sigma_s_per_m)
        else:
            z_s = complex(surface_impedance_fn(omega_lin))
        A, B, n_bc = _cavity_forms(
            V, 0.0, surface_impedance_ohm=z_s, omega_rad_s=omega_lin
        )

    if return_mode:
        # `EX-60`: additive — the eigenvalues are read exactly as without the flag.
        values, n_converged, vectors = _solve_pencil_nonhermitian(
            A, B, k0_sq, nev=nev, comm=comm, return_vectors=True
        )
    else:
        values, n_converged = _solve_pencil_nonhermitian(A, B, k0_sq, nev=nev, comm=comm)
    A.destroy()
    B.destroy()
    if n_converged == 0:
        raise RuntimeError("no eigenpair converged near the target")
    lam = complex(values[0])
    omega = C_0 * np.sqrt(lam)
    q = float(omega.real / (2.0 * abs(omega.imag))) if omega.imag != 0 else float("inf")

    mode_function = None
    if return_mode:
        # The eigenfunction of ``values[0]``; the owned block is filled and the
        # ghosts by scatter_forward (the ``solve_pec_cavity_modes`` pattern).
        mode_function = fem.Function(V, name="E_mode")
        local = mode_function.x.index_map.size_local * mode_function.x.block_size
        mode_function.x.array[:local] = vectors[0].getArray(readonly=True)
        mode_function.x.scatter_forward()
        for vec in vectors:
            vec.destroy()

    return ImpedanceWallMode(
        mode_function=mode_function,
        mesh=msh if return_mode else None,
        eigenvalue=lam,
        omega_rad_s=complex(omega),
        q=q,
        omega0_rad_s=omega0,
        n_converged=int(n_converged),
        n_cells=int(msh.topology.index_map(msh.topology.dim).size_global),
        n_dofs=int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs),
        n_constrained_dofs_local=int(n_bc),
        surface_impedance_ohm=z_s,
    )


def smallest_cavity_eigenvalues(
    edges: Sequence[float] = (1.0, 0.8, 0.6),
    divisions: Sequence[int] = (6, 5, 4),
    degree: int = 1,
    n_values: int = 8,
    comm: MPI.Comm = None,
) -> Tuple[np.ndarray, float]:
    """Return the ``n_values`` eigenvalues nearest zero, and the analytic ``k₁²``.

    This is the null-space diagnostic: the gradient modes must come back as a
    cluster at machine zero, cleanly separated from the first physical mode. A
    curl-curl assembly with a sign error, or a mass matrix that is not the plain
    ``L²`` form, generally contaminates this cluster.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    a, b, d = edges
    nx, ny, nz = divisions

    msh = dmesh.create_box(
        comm,
        [np.array([0.0, 0.0, 0.0]), np.array([a, b, d])],
        [nx, ny, nz],
        cell_type=dmesh.CellType.tetrahedron,
    )
    V = fem.functionspace(msh, ("N1curl", degree))

    lam_1 = float((2.0 * np.pi * analytic_cavity_frequencies(edges, 1)[0] / C_0) ** 2)
    bc_diagonal = 1.0e4 * lam_1

    A, B, _ = _cavity_forms(V, bc_diagonal)
    # Target just above zero: the gradient cluster is the nearest spectrum.
    values, _ = _solve_pencil(A, B, 1e-6 * lam_1, nev=n_values, comm=comm)
    A.destroy()
    B.destroy()
    return values[:n_values], lam_1


# ---------------------------------------------------------------------------
# `TH-17`: a general-mesh pencil beside ``_cavity_forms``.
#
# Purely additive.  ``_cavity_forms`` and every function above are untouched, so
# the box path (`TH-9`, `TH-14` step 1) is byte-identical.
#
# The problem solved here is the source-free, loaded, capacitively-loaded
# cavity on an arbitrary tetrahedral mesh:
#
#     ∫ (∇×u)·(∇×v) dx = k₀² [ ∫ ε_r,c u·v dx + Σ_s B_sheet,s(u, v) ]     (G1)
#
# with ``k₀² = ω²/c²``.  Two pieces are *imported*, never restated:
#
# * ``ε_r,c`` — the caller's complex relative permittivity per cell tag; a
#   conducting region is linearised at one ``ω_lin`` exactly as `TH-14` step 1
#   linearised the Leontovich wall (``ε_r − jσ/(ω_lin ε₀)`` on the ``e^{+jωt}``
#   convention this solver uses throughout), and the fixed point is a second
#   call at the found ``Re ω``.
# * ``B_sheet`` — the lumped-sheet term (L1) of
#   :func:`fem_em_solver.ports.lumped.lumped_port_bilinear_term`, *not* a
#   re-derivation.  (L1) is ``jωμ₀(1/R)∫(n×u)·(n×v) dS`` with
#   ``R = Z_p·w/h``; for a pure capacitor ``Z_p = 1/(jωC)`` this is
#   ``−ω²μ₀ C (h/w) ∫(n×u)·(n×v) dS`` — exactly ω², so
#
#       B_sheet = −(c²/ω_ref²) · L1(ω_ref)
#
#   is independent of ``ω_ref`` and moves the whole capacitor term onto the
#   right of (G1) as a surface *mass* term.  Building it through (L1) keeps the
#   restriction, the facet measure, the ``ufl.inner`` conjugation and
#   `PORT-14` step 3's κ width correction identical to the driven path.
#
# PEC is the ``_cavity_forms`` treatment: every exterior facet pinned by the
# large-diagonal elimination, ``A`` at ``bc_diagonal`` and ``B`` at 1.0.  On the
# `TH-15` step 3 hole mesh the coil's cavity wall *is* exterior, so this pins
# the outer box and the coil together — the same ``pec_facet_tags=None``
# reasoning the hole fixture documents.
# ---------------------------------------------------------------------------


def complex_relative_permittivity(
    epsilon_r: float, sigma_s_per_m: float, omega_rad_s: float
) -> complex:
    """``ε_r − jσ/(ω ε₀)`` — the ``e^{+jωt}`` lossy permittivity at one ω."""

    return complex(epsilon_r) - 1.0j * float(sigma_s_per_m) / (
        float(omega_rad_s) * EPSILON_0
    )


def general_mesh_pencil(
    V,
    epsilon_r_field,
    *,
    bc_diagonal: float,
    sheet_mass_forms: Sequence = (),
):
    """Assemble (G1)'s pencil ``(A, B)`` on a general mesh; see the note above.

    ``epsilon_r_field`` is any UFL-valid coefficient (a DG0 ``fem.Function`` is
    the intended shape); ``sheet_mass_forms`` are UFL forms already reduced to
    the ``k₀²`` side, i.e. ``−(c²/ω_ref²)·L1(ω_ref)`` per sheet.

    Returns ``(A, B, n_constrained_dofs_local)``.
    """

    if not np.issubdtype(PETSc.ScalarType, np.complexfloating):
        raise RuntimeError(
            "the general-mesh pencil needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode)"
        )

    msh = V.mesh
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    stiffness_ufl = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx
    mass_ufl = ufl.inner(epsilon_r_field * u, v) * ufl.dx
    for form in sheet_mass_forms:
        mass_ufl = mass_ufl + form

    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    boundary_facets = dmesh.exterior_facet_indices(msh.topology)
    boundary_dofs = fem.locate_dofs_topological(V, tdim - 1, boundary_facets)
    zero = fem.Function(V)
    zero.x.array[:] = 0.0
    bc = fem.dirichletbc(zero, boundary_dofs)

    A = assemble_matrix(fem.form(stiffness_ufl), bcs=[bc], diag=bc_diagonal)
    A.assemble()
    B = assemble_matrix(fem.form(mass_ufl), bcs=[bc], diag=1.0)
    B.assemble()
    return A, B, boundary_dofs.size


def solve_general_mesh_modes(
    V,
    epsilon_r_field,
    *,
    target_k0_sq: float,
    sheet_mass_forms: Sequence = (),
    nev: int = 6,
    bc_diagonal_factor: float = 1.0e4,
    comm: MPI.Comm = None,
    return_vectors: bool = False,
) -> Tuple[np.ndarray, int, int]:
    """Complex ``k₀²`` nearest ``target_k0_sq`` of (G1); GNHEP, shift-invert.

    Returns ``(eigenvalues ordered by |λ − target|, n_converged,
    n_constrained_dofs_local)``.  ``ω = c√λ`` and ``Q = Re ω/(2|Im ω|)``, the
    `TH-14` step 1 readout.  With ``return_vectors`` the eigenvectors follow
    as a fourth element, in the same order (`EX-76`; additive).
    """

    comm = V.mesh.comm if comm is None else comm
    A, B, n_bc = general_mesh_pencil(
        V,
        epsilon_r_field,
        bc_diagonal=float(bc_diagonal_factor) * float(target_k0_sq),
        sheet_mass_forms=sheet_mass_forms,
    )
    if return_vectors:
        values, n_converged, vectors = _solve_pencil_nonhermitian(
            A, B, float(target_k0_sq), nev=nev, comm=comm, return_vectors=True
        )
        A.destroy()
        B.destroy()
        return values, int(n_converged), int(n_bc), vectors
    values, n_converged = _solve_pencil_nonhermitian(
        A, B, float(target_k0_sq), nev=nev, comm=comm
    )
    A.destroy()
    B.destroy()
    return values, int(n_converged), int(n_bc)


# ---------------------------------------------------------------------------
# `GEO-33` step 0b: the Leontovich pencil on a supplied mesh.
#
# Purely additive: ``_cavity_forms`` and ``solve_impedance_wall_cavity_mode``
# are untouched.  On a mesh whose y = b/2 plane is duplicated (two lips, each
# an exterior facet with its own outward normal) the ordinary ``ds`` term of
# Jin §5.8.3 (impedance condition, Jin §1.5.3 eqs 1.54–1.57) already acts on
# both sides of the sheet.  ``sheet`` adds the single-valued interior form
#     γ Σ± ∫ (n±×u±)·(n±×v±) dS(tag),
# which on a tangentially continuous N1curl field is 2γ∫(n×u)·(n×v): the
# resistive sheet R = Z_s/2 of Jin §1.5.4 (eqs 1.60–1.63) — the control.
# ---------------------------------------------------------------------------


@dataclass
class MeshImpedanceModes:
    """Complex eigenvalues of the lossy pencil on a supplied mesh."""

    eigenvalues: np.ndarray
    n_converged: int
    n_dofs: int
    surface_impedance_ohm: complex
    #: Populated only with ``return_vectors=True`` (`EX-82`): one ``fem.Function``
    #: per converged eigenvector, in the |λ − k0²| order of ``eigenvalues``.
    mode_functions: list = None

    def q(self, i: int = 0) -> float:
        omega = C_0 * np.sqrt(complex(self.eigenvalues[i]))
        return float(omega.real / (2.0 * abs(omega.imag))) if omega.imag != 0 else float("inf")

    def omega(self, i: int = 0) -> complex:
        return complex(C_0 * np.sqrt(complex(self.eigenvalues[i])))


def solve_impedance_mesh_mode(
    msh,
    *,
    sigma_s_per_m: float,
    k0_sq: float,
    sheet=None,
    degree: int = 2,
    nev: int = 6,
    return_vectors: bool = False,
) -> MeshImpedanceModes:
    """Leontovich walls on every exterior facet of ``msh`` (± an interior sheet).

    ``Z_s = (1 + j)R_s(ω₀)`` with ``ω₀ = c√k0_sq`` (the linearisation, as in
    `TH-14`). ``sheet = (facet_tags, tag)`` adds the single-valued interior
    term over ``dS(tag)``. Eigenvalues are returned ordered by |λ − k0_sq|.
    """

    if not np.issubdtype(PETSc.ScalarType, np.complexfloating):
        raise RuntimeError(
            "a surface-impedance cavity needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode)"
        )
    V = fem.functionspace(msh, ("N1curl", degree))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    omega0 = float(np.sqrt(k0_sq) * C_0)
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega0, sigma_s_per_m)
    gamma = fem.Constant(msh, PETSc.ScalarType(1j * omega0 * MU_0 / z_s))
    n = ufl.FacetNormal(msh)
    a_form = (
        ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx
        + gamma * ufl.inner(ufl.cross(n, u), ufl.cross(n, v)) * ufl.ds
    )
    if sheet is not None:
        tags, tag = sheet
        dS = ufl.Measure("dS", domain=msh, subdomain_data=tags)
        for s in ("+", "-"):
            a_form = a_form + gamma * ufl.inner(
                ufl.cross(n(s), u(s)), ufl.cross(n(s), v(s))
            ) * dS(tag)
    A = assemble_matrix(fem.form(a_form))
    A.assemble()
    B = assemble_matrix(fem.form(ufl.inner(u, v) * ufl.dx))
    B.assemble()
    vectors = None
    if return_vectors:
        # `EX-82`: additive — eigenvalues are read exactly as without the flag.
        values, n_converged, vectors = _solve_pencil_nonhermitian(
            A, B, k0_sq, nev=nev, comm=msh.comm, return_vectors=True
        )
    else:
        values, n_converged = _solve_pencil_nonhermitian(A, B, k0_sq, nev=nev, comm=msh.comm)
    A.destroy()
    B.destroy()
    if n_converged == 0:
        raise RuntimeError("no eigenpair converged near the target")
    mode_functions = None
    if return_vectors:
        mode_functions = []
        for i, vec in enumerate(vectors):
            f = fem.Function(V, name=f"E_mode_{i}")
            local = f.x.index_map.size_local * f.x.block_size
            f.x.array[:local] = vec.getArray(readonly=True)
            f.x.scatter_forward()
            mode_functions.append(f)
            vec.destroy()
    return MeshImpedanceModes(
        eigenvalues=np.asarray(values, dtype=complex),
        n_converged=int(n_converged),
        n_dofs=int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs),
        surface_impedance_ohm=complex(z_s),
        mode_functions=mode_functions,
    )
