"""Validation test: layered-impedance reflection on a driven TEM line (TH-18 step 2).

Normal incidence on a planar stack — a volumetric dielectric slab (εr = 4,
0.2 m) backed by a conductor foil that is itself backed by free space —
against the closed-form transmission-line cascade (Jin §3.4.2 (3.87)–(3.89),
surface impedance Jin §1.5.3 (1.54)–(1.56)):

    Γ = (Z_in − η₀)/(Z_in + η₀),
    Z_in = layered_surface_impedance(ω, [DielectricLayer(4, 0.2), ConductorLayer(σ, t)], η₀).

In the FEM the slab is volumetric and **only the conductor** is the surface
term, Z_s = layered_surface_impedance(ω, [ConductorLayer(σ, t)], η₀), added
through ``TimeHarmonicSolver.solve(extra_bilinear_terms=…)`` as Jin §5.8.3's
third-kind term jωμ₀/Z_s ∫(n × E)·conj(n × W) dS (shaped like
``test_th14_birdcage_copper._leontovich_term``).

Fixture (PROJECT_PLAN §9 item 87): box 1.2 × 0.04 × 0.04 m, ``create_box``
(60, 2, 2), N1curl degree 2. Air x ∈ [0, 1.0], slab x ≥ 1.0. The drive face
x = 0 is pinned to E = ẑ, the z faces are PEC (TEM parallel-plate line, E_z),
the y faces are natural (PMC); the termination x = 1.2 carries the surface
term (pinned to 0 only in ``pec``). **Not `TH-6`'s fixture** — it pins the
exact total field on every face, so Γ from it would be circular.

E_z sampled at 41 off-lattice points x ∈ [0.15, 0.85]; LSQ
E_z = A e^{−jk₀x} + B e^{+jk₀x}; Γ_FEM = (B/A)e^{2jk₀·1.0} (e^{+jωt}).

Asserted: numpy anchors (n7)–(n9); premise (P) on ``leo`` ≤ 1e-3 in |Γ| and
∠Γ; (i) ``t5``, (ii) ``thin``, (iii) ``pec`` within the band frozen from W1;
negative control (rule (e)) bulk Leontovich vs the thin-layer closed form.
Printed only: copper 35 µm (PEC-limit reading), ``pmc``, slab-dropped ``thin``.
"""

import time

import dolfinx
import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem
from fem_em_solver.core.layered_impedance import (
    ConductorLayer,
    DielectricLayer,
    conductor_skin_depth_m,
    conductor_surface_resistance_ohm,
    layered_surface_impedance,
)
from fem_em_solver.core.time_harmonic import TimeHarmonicBoundaryCondition, TimeHarmonicSolver
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.constants import C_0, ETA_0, MU_0

pytestmark = pytest.mark.skipif(
    not np.issubdtype(PETSc.ScalarType, np.complexfloating),
    reason="driven TEM line needs the complex DolfinX build",
)

FREQ_HZ = 127.74e6
OMEGA = 2.0 * np.pi * FREQ_HZ
K0 = OMEGA / C_0
SIGMA = 10.0
EPS_SLAB = 4.0
L_AIR = 1.0
D_SLAB = 0.2
LENGTH = L_AIR + D_SLAB
WIDTH = 0.04
DIVISIONS = (60, 2, 2)
DEGREE = 2

AIR_TAG, SLAB_TAG = 1, 2
DRIVE_TAG, TERM_TAG, ZFACE_TAG = 11, 12, 13

N_SAMPLES = 41
X_SAMPLES = (0.15, 0.85)
YZ_FRACTION = 0.5137

COPPER_SIGMA = 5.8e7
COPPER_T = 35e-6

# Hand values (§9 item 87, review's two routes): (|Γ|, ∠Γ rad).
HAND_LEO = (0.9100, 1.5698)
HAND_THIN = (0.8382, 1.6375)
HAND_PEC_ARG = 1.6590
HAND_TOLERANCE = 1e-3  # (n8)
PEC_UNITY_TOLERANCE = 1e-14  # (n7)
T5_LEO_TOLERANCE = 1e-4  # (n9)
PREMISE_CAP = 1e-3  # (P), a-priori from TH-6's 0.059 % degree-1 β error
NEGATIVE_CONTROL_HAND = 0.093

# Band frozen from W1 (item 87): B = max(3·|Δ|_P, 1e-6), rounded up to one
# significant figure. Measured (P) on `leo`, 20260930T021219Z_TH-18-step2-W1.log:89:
#   d|G| = 7.185e-07, d_arg = 1.354e-05 rad  ->  3·d = 2.16e-06, 4.06e-05;
# band candidate printed at :91. Frozen after W1; never moved.
B_MAG = 3e-06  # 20260930T021219Z_TH-18-step2-W1.log:91
B_ARG = 5e-05  # rad, 20260930T021219Z_TH-18-step2-W1.log:91


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ceil_1sf(x):
    e = np.floor(np.log10(x))
    return float(np.ceil(x / 10.0**e - 1e-9) * 10.0**e)


def _delta():
    return conductor_skin_depth_m(OMEGA, SIGMA)


def _gamma_of(z_in):
    return complex((z_in - ETA_0) / (z_in + ETA_0))


def _slab():
    return DielectricLayer(EPS_SLAB, D_SLAB)


def _leo_zs():
    return (1.0 + 1.0j) * conductor_surface_resistance_ohm(OMEGA, SIGMA)


def _foil_zs(thickness, sigma=SIGMA):
    return layered_surface_impedance(OMEGA, [ConductorLayer(sigma, thickness)], ETA_0)


def closed_form_gamma(case):
    """Closed-form Γ at the slab front (x = 1.0) for each case."""
    delta = _delta()
    if case == "leo":
        z_in = layered_surface_impedance(OMEGA, [_slab()], _leo_zs())
    elif case == "t5":
        z_in = layered_surface_impedance(
            OMEGA, [_slab(), ConductorLayer(SIGMA, 5.0 * delta)], ETA_0
        )
    elif case == "thin":
        z_in = layered_surface_impedance(
            OMEGA, [_slab(), ConductorLayer(SIGMA, 0.5 * delta)], ETA_0
        )
    elif case == "pec":
        z_in = layered_surface_impedance(OMEGA, [_slab()], 0.0)
    elif case == "copper":
        z_in = layered_surface_impedance(
            OMEGA, [_slab(), ConductorLayer(COPPER_SIGMA, COPPER_T)], ETA_0
        )
    elif case == "pmc":
        # PMC backing: Z_in = Z_c / tanh(γd) (open-circuit load).
        z_c = ETA_0 / np.sqrt(EPS_SLAB)
        z_in = z_c / np.tanh(1.0j * K0 * np.sqrt(EPS_SLAB) * D_SLAB)
    elif case == "thin_noslab":
        z_in = _foil_zs(0.5 * _delta())
    else:
        raise KeyError(case)
    return _gamma_of(complex(z_in))


def _mag_arg(g):
    return abs(g), float(np.angle(g))


# ---------------------------------------------------------------------------
# numpy anchors (no FEM)
# ---------------------------------------------------------------------------


def test_n7_pec_unit_modulus():
    g = closed_form_gamma("pec")
    _say(f"[TH-18s2] n7 pec cf Gamma={g:.12f} |G|={abs(g):.16f} arg={np.angle(g):.6f}")
    assert abs(abs(g) - 1.0) <= PEC_UNITY_TOLERANCE
    assert abs(np.angle(g) - HAND_PEC_ARG) <= HAND_TOLERANCE


def test_n8_hand_values():
    delta = _delta()
    zc = _leo_zs()
    z_thin = _foil_zs(0.5 * delta)
    _say(
        f"[TH-18s2] n8 delta={delta * 1e3:.6f} mm Z_c={zc:.6f} Z_thin={z_thin:.6f} "
        f"sigma/(w eps0)={SIGMA / (OMEGA * 8.8541878128e-12):.2f}"
    )
    for case, (mag, arg) in (("leo", HAND_LEO), ("thin", HAND_THIN)):
        m, a = _mag_arg(closed_form_gamma(case))
        _say(f"[TH-18s2] n8 {case} cf |G|={m:.6f} arg={a:.6f} (hand {mag} / {arg})")
        assert abs(m - mag) <= HAND_TOLERANCE
        assert abs(a - arg) <= HAND_TOLERANCE


def test_n9_t5_equals_leontovich():
    d = abs(closed_form_gamma("t5") - closed_form_gamma("leo"))
    _say(f"[TH-18s2] n9 |G_t5 - G_leo| cf = {d:.3e}")
    assert d <= T5_LEO_TOLERANCE


# ---------------------------------------------------------------------------
# FEM
# ---------------------------------------------------------------------------


def _impedance_term(msh, facet_tags, z_s):
    """Jin §5.8.3's third-kind term on the termination — SpatialCoordinate-free."""
    ds_t = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags, subdomain_id=TERM_TAG)
    n = ufl.FacetNormal(msh)
    gamma = fem.Constant(msh, default_scalar_type(1j * OMEGA * MU_0 / complex(z_s)))

    def term(trial, test):
        return gamma * ufl.inner(ufl.cross(n, trial), ufl.cross(n, test)) * ds_t

    return term


def _meshtags(msh, dim, pairs):
    idx = np.concatenate([np.asarray(i, dtype=np.int32) for i, _ in pairs])
    val = np.concatenate([np.full(len(i), t, dtype=np.int32) for i, t in pairs])
    order = np.argsort(idx)
    idx, val = idx[order], val[order]
    assert np.unique(idx).size == idx.size
    return dolfinx.mesh.meshtags(msh, dim, idx, val)


_FIXTURE = {}


def _fixture():
    if _FIXTURE:
        return _FIXTURE
    comm = MPI.COMM_WORLD
    msh = dolfinx.mesh.create_box(
        comm, [np.array([0.0, 0.0, 0.0]), np.array([LENGTH, WIDTH, WIDTH])], list(DIVISIONS)
    )
    tdim = msh.topology.dim
    fdim = tdim - 1
    slab = dolfinx.mesh.locate_entities(msh, tdim, lambda x: x[0] >= L_AIR - 1e-9)
    n_local = msh.topology.index_map(tdim).size_local
    all_cells = np.arange(n_local, dtype=np.int32)
    slab = slab[slab < n_local]
    air = np.setdiff1d(all_cells, slab)
    cell_tags = _meshtags(msh, tdim, [(air, AIR_TAG), (slab, SLAB_TAG)])

    msh.topology.create_connectivity(fdim, tdim)
    drive = dolfinx.mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], 0.0))
    term = dolfinx.mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], LENGTH))
    zf = dolfinx.mesh.locate_entities_boundary(
        msh, fdim, lambda x: np.isclose(x[2], 0.0) | np.isclose(x[2], WIDTH)
    )
    zf = np.setdiff1d(zf, np.union1d(drive, term))
    facet_tags = _meshtags(msh, fdim, [(drive, DRIVE_TAG), (term, TERM_TAG), (zf, ZFACE_TAG)])

    n_cells = comm.allreduce(n_local, op=MPI.SUM)
    n_slab = comm.allreduce(int(slab.size), op=MPI.SUM)
    n_term = comm.allreduce(
        int(np.count_nonzero(np.asarray(facet_tags.values) == TERM_TAG)), op=MPI.SUM
    )
    _say(f"[TH-18s2] mesh cells={n_cells} slab_cells={n_slab} term_facets(rank-summed)={n_term}")
    assert n_slab == n_cells // 6  # 0.2 of 1.2 m
    xs = np.linspace(*X_SAMPLES, N_SAMPLES)
    pts = np.column_stack(
        [xs, np.full_like(xs, YZ_FRACTION * WIDTH), np.full_like(xs, YZ_FRACTION * WIDTH)]
    )
    _FIXTURE.update(
        msh=msh, cell_tags=cell_tags, facet_tags=facet_tags, xs=xs, pts=pts, n_cells=n_cells
    )
    return _FIXTURE


def _drive(x):
    vals = np.zeros((3, x.shape[1]), dtype=np.complex128)
    vals[2] = (x[0] < 0.6).astype(np.float64)
    return vals


_SOLVES = {}


def fem_gamma(case, return_fields=False):
    """Solve one termination and return Γ_FEM referenced at x = 1.0.

    ``return_fields=True`` (EX-70, additive): skip the cache read, solve, store Γ in
    ``_SOLVES`` as always and return ``(Γ, fields)``.  The default path is unchanged.
    """
    if case in _SOLVES and not return_fields:
        return _SOLVES[case]
    fx = _fixture()
    msh, facet_tags = fx["msh"], fx["facet_tags"]
    pec_tags = [DRIVE_TAG, ZFACE_TAG] + ([TERM_TAG] if case == "pec" else [])
    delta = _delta()
    zs = {
        "leo": _leo_zs(),
        "t5": _foil_zs(5.0 * delta) if case == "t5" else None,
        "thin": _foil_zs(0.5 * delta) if case == "thin" else None,
        "copper": _foil_zs(COPPER_T, COPPER_SIGMA) if case == "copper" else None,
        "pec": None,
        "pmc": None,
    }[case]
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQ_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
        material_map={SLAB_TAG: HomogeneousMaterial(sigma=0.0, epsilon_r=EPS_SLAB)},
        cell_tags=fx["cell_tags"],
        facet_tags=facet_tags,
        boundary_condition=TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A,
        dirichlet_e_field=_drive,
        pec_facet_tags=pec_tags,
    )
    t = time.perf_counter()
    solver = TimeHarmonicSolver(problem, degree=DEGREE)
    terms = [] if zs is None else [_impedance_term(msh, facet_tags, zs)]
    fields = solver.solve(extra_bilinear_terms=terms or None)
    vals, mask = evaluate_vector_field_parallel(fields.e_complex, fx["pts"])
    assert bool(np.all(mask))
    ez = vals[:, 2]
    xs = fx["xs"]
    m = np.column_stack([np.exp(-1j * K0 * xs), np.exp(1j * K0 * xs)])
    (a, b), *_ = np.linalg.lstsq(m, ez, rcond=None)
    resid = float(np.linalg.norm(m @ np.array([a, b]) - ez) / np.linalg.norm(ez))
    g = complex(b / a * np.exp(2j * K0 * L_AIR))
    ndofs = (
        solver.function_space().dofmap.index_map.size_global
        * solver.function_space().dofmap.index_map_bs
    )
    _say(
        f"[TH-18s2] FEM {case:6s} Z_s={zs} dofs={ndofs} A={a:.6f} B={b:.6f} "
        f"|1+G(0)|={abs(1 + b / a):.4f} lsq_resid={resid:.2e} Gamma={g:.8f} "
        f"|G|={abs(g):.8f} arg={np.angle(g):.8f} ({time.perf_counter() - t:.1f}s)"
    )
    _SOLVES[case] = g
    if return_fields:
        return g, fields
    return g


def _deltas(case):
    g_fem = fem_gamma(case)
    g_cf = closed_form_gamma(case)
    d_mag = abs(abs(g_fem) - abs(g_cf))
    d_arg = abs(float(np.angle(g_fem / g_cf)))
    _say(
        f"[TH-18s2] {case:6s} cf |G|={abs(g_cf):.8f} arg={np.angle(g_cf):.8f} "
        f"| d|G|={d_mag:.3e} d_arg={d_arg:.3e} rad"
    )
    return d_mag, d_arg


def test_premise_leontovich_reflection():
    d_mag, d_arg = _deltas("leo")
    nc = abs(fem_gamma("leo") - closed_form_gamma("thin"))
    _say(f"[TH-18s2] premise |G_leo,FEM - G_thin,cf| = {nc:.6f} (hand {NEGATIVE_CONTROL_HAND})")
    b_mag = _ceil_1sf(max(3.0 * d_mag, 1e-6))
    b_arg = _ceil_1sf(max(3.0 * d_arg, 1e-6))
    _say(f"[TH-18s2] band candidate B_mag={b_mag:.0e} B_arg={b_arg:.0e} rad")
    assert d_mag <= PREMISE_CAP
    assert d_arg <= PREMISE_CAP


def _bands():
    if B_MAG is None or B_ARG is None:
        pytest.fail("band not frozen — W1 must run first and its values be written into the module")
    return B_MAG, B_ARG


@pytest.mark.parametrize("case", ["t5", "thin", "pec"])
def test_gate_within_frozen_band(case):
    b_mag, b_arg = _bands()
    d_mag, d_arg = _deltas(case)
    assert d_mag <= b_mag
    assert d_arg <= b_arg


def test_negative_control_bulk_leontovich_vs_thin():
    b_mag, b_arg = _bands()
    nc = abs(fem_gamma("leo") - closed_form_gamma("thin"))
    b = max(b_mag, b_arg)
    _say(
        f"[TH-18s2] negative control |G_leo,FEM - G_thin,cf| = {nc:.6f}; floor 10*B = {10 * b:.1e}"
    )
    assert nc > 10.0 * b
    assert abs(nc - NEGATIVE_CONTROL_HAND) <= b  # item 87: ceiling fixed at 0.093 ± B


def test_printed_only_readings():
    g_pec = fem_gamma("pec")
    g_cu = fem_gamma("copper")
    _say(
        f"[TH-18s2] printed copper: |G_cu,FEM - G_pec,FEM| = {abs(g_cu - g_pec):.3e}; "
        f"cf |G_cu - G_pec| = {abs(closed_form_gamma('copper') - closed_form_gamma('pec')):.3e}"
    )
    g_pmc = fem_gamma("pmc")
    g_pmc_cf = closed_form_gamma("pmc")
    _say(
        f"[TH-18s2] printed pmc: FEM |G|={abs(g_pmc):.6f} arg={np.angle(g_pmc):.6f}; "
        f"cf |G|={abs(g_pmc_cf):.6f} arg={np.angle(g_pmc_cf):.6f}"
    )
    g_ns = closed_form_gamma("thin_noslab")
    _say(f"[TH-18s2] printed slab-dropped thin cf |G|={abs(g_ns):.6f} arg={np.angle(g_ns):.6f}")
