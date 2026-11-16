"""`TH-5` step 2b (§9 item 290): the mixed wall on `PORT-1`'s two-torus box.

Two solves at `PORT-1`'s smallest 3b-xi padding W = 0.08 m on one mesh (PEC, the
premise; natural/magnetic wall, the measurement) and two printed natural rungs
at 0.10 / 0.12 m.  The mixed impedance is the post-solve combination
Z^mix = ⅓ Z^PEC + ⅔ Z^NAT (the §7 entry's "mixed average"); it is not a
boundary condition and not a face split.

Weights from the sphere model (stated before the run): a centred magnetic
dipole in a sphere of radius R sees a uniform reaction field of −2m/(4πR³)
under PEC and +m/(4πR³) under a magnetic wall, so B^NAT = −½ B^PEC on the
l = 1 term and w_PEC = ⅓, w_NAT = ⅔ cancel it (2w_PEC − w_NAT = 0,
w_PEC + w_NAT = 1).

Box term B^X(W) = d^X(W) − D∞, d^X = 100·(Im Z₂₁^X / ωM₁₂ − 1) pp, D∞ = D∞^PEC
from `_free_exponent_fit` over `PORT1_PADDING_RUNGS_M` /
`PORT1_PADDING_DEFICITS_PP` (`GEO-4` (A0)).  Weights sum to 1, so
B^mix = ⅓ B^PEC + ⅔ B^NAT = d^mix − D∞.

Asserted: (P0) PEC 0.08 m reproduces `P0_RUNG_PP` within `P0_BAND_PP`;
(E) |Re Z₂₁| ≤ IDENTITY_TOLERANCE·|Z₂₁| under both walls; (S) B^NAT(0.08) > 0
against B^PEC(0.08) < 0, and (c1) d^NAT, d^PEC differ by > P0_BAND_PP;
(M) |B^mix(0.08)| ≤ |B^PEC(0.08)| / MIX_BAR_DIVISOR; (c2) swapped weights
⅔ PEC + ⅓ NAT miss the bar.  Everything else printed only.

The solve is a copy of `_solve_reaction_z`'s mesh / current / solve / reaction
lines (`test_port_reaction_impedance.py`), restricted to wall ∈ {"pec",
"natural"}; that gated module is not edited.  (P0) proves the copy is the
fixture.

Brief: docs/planning/chunks/TH-5.md, "§9 item 290".
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import (
    PORT1_PADDING_DEFICITS_PP,
    PORT1_PADDING_RUNGS_M,
    MeshGenerator,
    _free_exponent_fit,
)

from tests.complex_mode import complex_only
from tests.validation.test_geo4_box_padding_heuristic import P0_BAND_PP, P0_RUNG_PP
from tests.validation.test_port_reaction_impedance import (
    CURRENT_A,
    FREQUENCY_HZ,
    H_FAR,
    H_WIRE,
    IDENTITY_TOLERANCE,
    MAJOR_RADIUS,
    MINOR_RADIUS,
    OMEGA,
    SEPARATION,
    TORUS_TAGS,
    _azimuthal_current_density,
    _reaction,
    _tag_volume,
    mutual_inductance,
)

# The weekly's weights and divisor (sphere model, stated before the run; never fitted).
W_PEC = 1.0 / 3.0
W_NAT = 2.0 / 3.0
MIX_BAR_DIVISOR = 5.0

PREMISE_PADDING = 0.08
PRINTED_NAT_PADDINGS = (0.10, 0.12)


def _root_print(msg: str) -> None:
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _d_inf_pec() -> float:
    d_inf, _p = _free_exponent_fit(PORT1_PADDING_RUNGS_M, PORT1_PADDING_DEFICITS_PP)
    return d_inf


def _b_pec_fit(w: float) -> float:
    i = list(PORT1_PADDING_RUNGS_M).index(w)
    return PORT1_PADDING_DEFICITS_PP[i] - _d_inf_pec()


def _build_mesh(padding: float):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = MeshGenerator.two_torus_domain(
        separation=SEPARATION,
        major_radius=MAJOR_RADIUS,
        minor_radius=MINOR_RADIUS,
        resolution=H_FAR,
        air_padding=padding,
        wire_resolution=H_WIRE,
        far_resolution=H_FAR,
        comm=comm,
    )
    t_mesh = time.perf_counter() - t0
    ncells = msh.topology.index_map(msh.topology.dim).size_global
    _root_print(
        f"[TH-5 2b] mesh padding {padding} m: {ncells} cells, {t_mesh:.1f} s"
    )
    return msh, cell_tags, ncells


def _solve_z_column(msh, cell_tags, ncells: int, padding: float, wall: str) -> dict:
    """Copy of `_solve_reaction_z`'s current / solve / reaction lines, tag 1 driven.

    ``wall`` ∈ {"pec", "natural"}; natural = no Dirichlet rows, no surface term.
    """
    assert wall in ("pec", "natural")
    comm = MPI.COMM_WORLD

    j_magnitude = CURRENT_A / (np.pi * MINOR_RADIUS**2)
    volumes = [_tag_volume(msh, cell_tags, t, comm) for t in TORUS_TAGS]
    currents = [j_magnitude * v / (2.0 * np.pi * MAJOR_RADIUS) for v in volumes]

    driven_tag = TORUS_TAGS[0]
    col = 0
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=cell_tags,
        material_map=None,
        boundary_condition="pec_zero_tangential_a" if wall == "pec" else "natural",
    )
    solver = TimeHarmonicSolver(problem, degree=1)
    comm.Barrier()
    t0 = time.perf_counter()
    fields = solver.solve(
        current_density=_azimuthal_current_density(j_magnitude),
        subdomain_ids=[driven_tag],
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0

    z = np.zeros(2, dtype=complex)
    for row, test_tag in enumerate(TORUS_TAGS):
        z[row] = _reaction(
            msh, cell_tags, fields.e_complex, test_tag, currents[col], currents[row], comm
        )

    omega_m12 = OMEGA * mutual_inductance(MAJOR_RADIUS, MAJOR_RADIUS, SEPARATION)
    deficit = 100.0 * (z[1].imag / omega_m12 - 1.0)
    # Lossless reaction/energy residual (printed): Re Z / |Z| on both entries.
    res11 = abs(z[0].real) / abs(z[0])
    res21 = abs(z[1].real) / abs(z[1])
    _root_print(
        f"[TH-5 2b] wall={wall} padding {padding} m: {ncells} cells, "
        f"solve {t_solve:.1f} s; Z11 = {z[0].real:+.6e}{z[0].imag:+.6e}j Ohm, "
        f"Z21 = {z[1].real:+.6e}{z[1].imag:+.6e}j Ohm, omega*M12 = {omega_m12:+.6e} Ohm, "
        f"deficit {deficit:+.4f} pp; |Re Z11|/|Z11| {res11:.3e}, "
        f"|Re Z21|/|Z21| {res21:.3e}; I1 = {currents[0]:.6f} A"
    )
    return {
        "z11": complex(z[0]),
        "z21": complex(z[1]),
        "omega_m12": omega_m12,
        "deficit_pp": deficit,
        "ncells": ncells,
        "t_solve": t_solve,
    }


@pytest.fixture(scope="module")
def runs():
    out = {}
    msh, cell_tags, ncells = _build_mesh(PREMISE_PADDING)
    out[("pec", PREMISE_PADDING)] = _solve_z_column(msh, cell_tags, ncells, PREMISE_PADDING, "pec")
    out[("natural", PREMISE_PADDING)] = _solve_z_column(
        msh, cell_tags, ncells, PREMISE_PADDING, "natural"
    )
    del msh, cell_tags
    for w in PRINTED_NAT_PADDINGS:
        msh, cell_tags, ncells = _build_mesh(w)
        out[("natural", w)] = _solve_z_column(msh, cell_tags, ncells, w, "natural")
        del msh, cell_tags
    return out


def _mix_deficit(pec: dict, nat: dict, w_pec: float, w_nat: float) -> float:
    """Deficit of Z^mix = w_pec Z^PEC + w_nat Z^NAT (same ωM₁₂, so linear in Z)."""
    z21 = w_pec * pec["z21"] + w_nat * nat["z21"]
    return 100.0 * (z21.imag / pec["omega_m12"] - 1.0)


def test_sphere_model_weights_numpy():
    """Printed check of the l-pole reaction coefficients: PEC (l+1)/l, magnetic −1.

    Asserted only as the algebraic identity the weights are built from
    (2·W_PEC − W_NAT = 0, W_PEC + W_NAT = 1) — exact.
    """
    for l in (1, 2, 3):
        c_pec, c_nat = (l + 1) / l, -1.0
        resid = W_PEC * c_pec + W_NAT * c_nat
        _root_print(
            f"[TH-5 2b sphere] l={l}: PEC coeff {c_pec:+.4f}, NAT coeff {c_nat:+.4f}; "
            f"weighted residual / PEC term = {resid / c_pec:+.4f} "
            f"(closed form (1-l)/(3l) = {(1 - l) / (3 * l):+.4f})"
        )
    assert abs(2.0 * W_PEC - W_NAT) <= 1e-15
    assert abs(W_PEC + W_NAT - 1.0) <= 1e-15


@complex_only
def test_p0_premise_pec_rung(runs):
    d = runs[("pec", PREMISE_PADDING)]["deficit_pp"]
    delta = abs(d - P0_RUNG_PP)
    _root_print(
        f"[TH-5 2b P0] PEC 0.08 m deficit {d:+.4f} pp vs {P0_RUNG_PP:+.4f} pp: "
        f"|delta| {delta:.4f} pp (band {P0_BAND_PP})"
    )
    assert delta <= P0_BAND_PP


@complex_only
def test_e_lossless_identity_both_walls(runs):
    rels = {}
    for wall in ("pec", "natural"):
        r = runs[(wall, PREMISE_PADDING)]
        rels[wall] = abs(r["z21"].real) / abs(r["z21"])
        _root_print(
            f"[TH-5 2b E] {wall} 0.08 m: |Re Z21|/|Z21| = {rels[wall]:.3e} "
            f"(tol {IDENTITY_TOLERANCE})"
        )
    for w in PRINTED_NAT_PADDINGS:
        r = runs[("natural", w)]
        _root_print(
            f"[TH-5 2b E] natural {w} m: |Re Z21|/|Z21| = "
            f"{abs(r['z21'].real) / abs(r['z21']):.3e} (printed)"
        )
    for wall in ("pec", "natural"):
        assert rels[wall] <= IDENTITY_TOLERANCE


@complex_only
def test_s_m_c2_mixed_wall(runs):
    d_inf = _d_inf_pec()
    pec = runs[("pec", PREMISE_PADDING)]
    nat = runs[("natural", PREMISE_PADDING)]
    d_pec, d_nat = pec["deficit_pp"], nat["deficit_pp"]
    b_pec = d_pec - d_inf
    b_nat = d_nat - d_inf
    r_nat = abs(b_nat) / abs(b_pec)
    d_mix = _mix_deficit(pec, nat, W_PEC, W_NAT)
    b_mix = d_mix - d_inf
    bar = abs(b_pec) / MIX_BAR_DIVISOR
    d_swap = _mix_deficit(pec, nat, W_NAT, W_PEC)
    b_swap = d_swap - d_inf
    d_half = _mix_deficit(pec, nat, 0.5, 0.5)
    b_half = d_half - d_inf
    w_star = r_nat / (1.0 + r_nat)
    literal_bar = abs(d_inf) / MIX_BAR_DIVISOR

    _root_print(
        f"[TH-5 2b S] W 0.08 m: D_inf^PEC {d_inf:+.10f} pp; d^PEC {d_pec:+.4f} pp, "
        f"d^NAT {d_nat:+.4f} pp; B^PEC {b_pec:+.4f} pp, B^NAT {b_nat:+.4f} pp; "
        f"r_NAT {r_nat:.4f}"
    )
    _root_print(
        f"[TH-5 2b c1] |d^NAT - d^PEC| = {abs(d_nat - d_pec):.4f} pp "
        f"(must exceed P0_BAND_PP {P0_BAND_PP}); bitwise equal: {d_nat == d_pec}"
    )
    _root_print(
        f"[TH-5 2b M] d^mix {d_mix:+.4f} pp, B^mix {b_mix:+.4f} pp; "
        f"bar |B^PEC|/{MIX_BAR_DIVISOR:g} = {bar:.4f} pp; pass: {abs(b_mix) <= bar}; "
        f"literal-reading bar |D_inf^PEC|/{MIX_BAR_DIVISOR:g} = {literal_bar:.4f} pp, "
        f"met: {abs(b_mix) <= literal_bar} (printed only)"
    )
    _root_print(
        f"[TH-5 2b c2] swapped (2/3 PEC + 1/3 NAT): B^swap {b_swap:+.4f} pp, "
        f"|B^swap|/bar = {abs(b_swap) / bar:.4f}; equal (1/2,1/2): B^half "
        f"{b_half:+.4f} pp, |B^half|/bar = {abs(b_half) / bar:.4f} (printed)"
    )
    _root_print(
        f"[TH-5 2b w*] zeroing weights on this cube: w_PEC* = r/(1+r) = {w_star:.4f} "
        f"(vs W_PEC {W_PEC:.4f}), w_NAT* = {1.0 - w_star:.4f} (printed only)"
    )
    # Printed rungs: natural at 0.10 / 0.12 m against the fit-referenced B^PEC.
    _root_print(
        f"[TH-5 2b rungs] W 0.08 m: r_NAT {r_nat:.4f}, B^mix {b_mix:+.4f} pp "
        f"(B^PEC from this window's solve {b_pec:+.4f}; fit-referenced "
        f"{_b_pec_fit(PREMISE_PADDING):+.4f})"
    )
    for w in PRINTED_NAT_PADDINGS:
        d_n = runs[("natural", w)]["deficit_pp"]
        b_n = d_n - d_inf
        b_p = _b_pec_fit(w)
        r_w = abs(b_n) / abs(b_p)
        b_m = W_PEC * b_p + W_NAT * b_n
        _root_print(
            f"[TH-5 2b rungs] W {w} m (fit-referenced B^PEC {b_p:+.4f} pp): "
            f"d^NAT {d_n:+.4f} pp, B^NAT {b_n:+.4f} pp, r_NAT {r_w:.4f}, "
            f"B^mix {b_m:+.4f} pp, |B^mix|/|B^PEC| {abs(b_m) / abs(b_p):.4f}, "
            f"w_PEC* {r_w / (1.0 + r_w):.4f}"
        )

    # (S) sign flip of the natural wall; (c1) wall actually switched.
    assert b_pec < 0.0
    assert b_nat > 0.0
    assert abs(d_nat - d_pec) > P0_BAND_PP
    # (M) the weekly's gate.
    assert abs(b_mix) <= bar
    # (c2) swapped weights miss the bar (conditional on (M) green, reached only then).
    assert abs(b_swap) > bar
