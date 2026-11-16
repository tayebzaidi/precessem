"""`TH-5` step 2a (§9 item 277): the first-order ABC on `PORT-1`'s two-torus box.

Three solves at `PORT-1`'s 3b-xi paddings W = 0.08 / 0.10 / 0.12 m with the PEC
wall replaced by the first-order ABC (Jin 3e §9.2 (9.44), s = 1), plus a PEC
premise solve and a natural-wall control at 0.08 m.  At 10 MHz the box sits at
k₀·half-width ≈ 0.022–0.035, where the ABC is a reflector equal to the natural
(magnetic) wall to O(k₀R).  The sphere model's static reaction fields (+2 under
PEC, −1 under a magnetic wall) predict the box term flips sign and halves.

Box term B(W) = deficit(W) − D∞^PEC, D∞^PEC from `_free_exponent_fit` over
`PORT1_PADDING_RUNGS_M` / `PORT1_PADDING_DEFICITS_PP` (`GEO-4` (A0)).

Asserted: (P0) the copied solver reproduces the PEC 0.08 m rung within
`P0_BAND_PP`; (S) B^ABC > 0 at all three rungs (B^PEC < 0); (E) Re Z₁₁ > 0
under ABC and |Re Z₂₁| ≤ IDENTITY_TOLERANCE·|Z₂₁| under PEC.  Controls (c1)
PEC in the ABC's place — via (S) and (P0); (c2) conjugated γ — via (E)'s
sign; (c3) natural wall at 0.08 m — printed only.  The entry's original ≥ 5×
bar is printed, never asserted (the banked negative).

Printed (record): (R) r(W) = |B^ABC|/|B^PEC| against the sphere model's 0.5 —
retired as an assertion by the 2026-10-07 weekly (a changed quantity, not a
loosened band; §9 item 320).

Brief: docs/planning/chunks/TH-5.md, "§9 item 277" and "§9 item 320".
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from dolfinx import fem

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.core.boundary_terms import first_order_abc_term
from fem_em_solver.io.mesh import (
    PORT1_PADDING_DEFICITS_PP,
    PORT1_PADDING_RUNGS_M,
    MeshGenerator,
    _free_exponent_fit,
)
from fem_em_solver.utils.constants import C_0, MU_0

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

# The entry's original "≥ 5×" bar — printed only (the banked negative).
ORIGINAL_FALL_FACTOR = 5.0
OUTER_BOUNDARY_TAG = 1  # `two_torus_domain` facet tag 1 (io/mesh.py:1353–1354)
ABC_PADDINGS = tuple(PORT1_PADDING_RUNGS_M)  # 0.08 / 0.10 / 0.12 m
PREMISE_PADDING = 0.08
ETA_0 = MU_0 * C_0


def _root_print(msg: str) -> None:
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _d_inf_pec() -> float:
    d_inf, _p = _free_exponent_fit(PORT1_PADDING_RUNGS_M, PORT1_PADDING_DEFICITS_PP)
    return d_inf


def _solve_z_column(padding: float, wall: str, return_fields: bool = False) -> dict:
    """Copy of `_solve_reaction_z`'s mesh / current / solve / reaction lines.

    One mesh at ``padding``, tag 1 driven, ``wall`` ∈ {"pec", "natural", "abc"}.
    "abc" = natural plus the first-order ABC term on the outer boundary.
    """
    assert wall in ("pec", "natural", "abc")
    comm = MPI.COMM_WORLD
    t_mesh = time.perf_counter()
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
    t_mesh = time.perf_counter() - t_mesh
    tdim = msh.topology.dim
    ncells = msh.topology.index_map(tdim).size_global

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
    extra = None
    if wall == "abc":
        extra = [first_order_abc_term(msh, facet_tags, OUTER_BOUNDARY_TAG, OMEGA)]
    solver = TimeHarmonicSolver(problem, degree=1)
    comm.Barrier()
    t0 = time.perf_counter()
    fields = solver.solve(
        current_density=_azimuthal_current_density(j_magnitude),
        subdomain_ids=[driven_tag],
        extra_bilinear_terms=extra,
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0

    z = np.zeros(2, dtype=complex)
    for row, test_tag in enumerate(TORUS_TAGS):
        z[row] = _reaction(
            msh, cell_tags, fields.e_complex, test_tag, currents[col], currents[row], comm
        )

    # Printed only: ∮_{tag 1} |n × E_h|² ds / η₀ (rank-reduced, real part after).
    n = ufl.FacetNormal(msh)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags)
    e = fields.e_complex
    nxe = ufl.cross(n, e)
    local = fem.assemble_scalar(fem.form(ufl.inner(nxe, nxe) * ds(OUTER_BOUNDARY_TAG)))
    boundary_power = float(np.real(comm.allreduce(local, op=MPI.SUM))) / ETA_0

    omega_m12 = OMEGA * mutual_inductance(MAJOR_RADIUS, MAJOR_RADIUS, SEPARATION)
    deficit = 100.0 * (z[1].imag / omega_m12 - 1.0)
    _root_print(
        f"[TH-5 2a] wall={wall} padding {padding} m: {ncells} cells, mesh {t_mesh:.1f} s, "
        f"solve {t_solve:.1f} s; Z11 = {z[0].real:+.6e}{z[0].imag:+.6e}j Ohm, "
        f"Z21 = {z[1].real:+.6e}{z[1].imag:+.6e}j Ohm, omega*M12 = {omega_m12:+.6e} Ohm, "
        f"deficit {deficit:+.4f} pp; I1 = {currents[0]:.6f} A, "
        f"oint|n x E|^2/eta0 = {boundary_power:.6e} W"
    )
    out = {
        "z11": complex(z[0]),
        "z21": complex(z[1]),
        "deficit_pp": deficit,
        "i1": currents[0],
        "boundary_power": boundary_power,
        "ncells": ncells,
    }
    if return_fields:
        # Additive (EX-103 rule (a) lift): the field Function only, never `fields`.
        out["omega_m12"] = omega_m12
        out["e"] = fields.e_complex
        out["msh"] = msh
        out["cell_tags"] = cell_tags
        out["facet_tags"] = facet_tags
    return out


@pytest.fixture(scope="module")
def runs():
    out = {("pec", PREMISE_PADDING): _solve_z_column(PREMISE_PADDING, "pec")}
    for w in ABC_PADDINGS:
        out[("abc", w)] = _solve_z_column(w, "abc")
    out[("natural", PREMISE_PADDING)] = _solve_z_column(PREMISE_PADDING, "natural")
    return out


def _b_pec(w: float) -> float:
    i = list(PORT1_PADDING_RUNGS_M).index(w)
    return PORT1_PADDING_DEFICITS_PP[i] - _d_inf_pec()


@complex_only
def test_p0_premise_pec_rung(runs):
    d = runs[("pec", PREMISE_PADDING)]["deficit_pp"]
    delta = abs(d - P0_RUNG_PP)
    _root_print(
        f"[TH-5 2a P0] PEC 0.08 m deficit {d:+.4f} pp vs {P0_RUNG_PP:+.4f} pp: "
        f"|delta| {delta:.4f} pp (band {P0_BAND_PP})"
    )
    assert delta <= P0_BAND_PP


@complex_only
def test_e_abc_absorbs_and_pec_structural_zero(runs):
    pec = runs[("pec", PREMISE_PADDING)]
    rel = abs(pec["z21"].real) / abs(pec["z21"])
    _root_print(f"[TH-5 2a E] PEC |Re Z21|/|Z21| = {rel:.3e} (tol {IDENTITY_TOLERANCE})")
    for w in ABC_PADDINGS:
        r = runs[("abc", w)]
        lhs = r["z11"].real * r["i1"] ** 2
        resid = abs(lhs - r["boundary_power"]) / abs(r["boundary_power"])
        _root_print(
            f"[TH-5 2a E] ABC {w} m: Re Z11 = {r['z11'].real:+.6e} Ohm; "
            f"Re Z11*I1^2 = {lhs:.6e} W vs oint|nxE|^2/eta0 = {r['boundary_power']:.6e} W, "
            f"rel residual {resid:.3e} (predicted <= 1e-6, printed only)"
        )
    nat = runs[("natural", PREMISE_PADDING)]
    _root_print(f"[TH-5 2a E] natural 0.08 m: Re Z11 = {nat['z11'].real:+.6e} Ohm (printed)")
    assert rel <= IDENTITY_TOLERANCE
    for w in ABC_PADDINGS:
        assert runs[("abc", w)]["z11"].real > 0.0


@complex_only
def test_s_r_box_term_sign_and_ratio(runs):
    d_inf = _d_inf_pec()
    ratios, b_abc_all = [], []
    for w in ABC_PADDINGS:
        b_pec = _b_pec(w)
        b_abc = runs[("abc", w)]["deficit_pp"] - d_inf
        r = abs(b_abc) / abs(b_pec)
        ratios.append(r)
        b_abc_all.append(b_abc)
        _root_print(
            f"[TH-5 2a S/R] W {w} m: deficit^ABC {runs[('abc', w)]['deficit_pp']:+.4f} pp, "
            f"B^ABC {b_abc:+.4f} pp, B^PEC {b_pec:+.4f} pp, r = {r:.4f} "
            f"(record; (R) retired, weekly 2026-10-07; sphere model 0.5); original bar |B^ABC| <= |B^PEC|/{ORIGINAL_FALL_FACTOR:g}: "
            f"{abs(b_abc) <= abs(b_pec) / ORIGINAL_FALL_FACTOR} (printed only)"
        )
    # (c1) separation of the ABC from a PEC-like broken run at 0.08 m (printed).
    _root_print(
        f"[TH-5 2a c1] |B^ABC - B^PEC| at 0.08 m = "
        f"{abs(b_abc_all[0] - _b_pec(PREMISE_PADDING)):.4f} pp"
    )
    # (c3) natural wall — printed only.
    nat = runs[("natural", PREMISE_PADDING)]["deficit_pp"]
    abc8 = runs[("abc", PREMISE_PADDING)]["deficit_pp"]
    _root_print(
        f"[TH-5 2a c3] |deficit^ABC - deficit^NAT| at 0.08 m = {abs(abc8 - nat):.6f} pp "
        f"(predicted <= {0.05 * abs(_b_pec(PREMISE_PADDING)):.4f} pp, printed only)"
    )
    # Free-exponent fit through the ABC ladder — printed only.
    abc_def = tuple(runs[("abc", w)]["deficit_pp"] for w in ABC_PADDINGS)
    d_abc, p_abc = _free_exponent_fit(ABC_PADDINGS, abc_def)
    if d_abc is None:
        _root_print(f"[TH-5 2a fit] free-exponent fit failed on ABC ladder {abc_def}")
    else:
        _root_print(
            f"[TH-5 2a fit] D_inf^ABC = {d_abc:+.6f} pp, p^ABC = {p_abc:.4f}, "
            f"|D_inf^ABC - D_inf^PEC| = {abs(d_abc - d_inf):.6f} pp (D_inf^PEC {d_inf:+.6f})"
        )
    for w in ABC_PADDINGS:
        assert _b_pec(w) < 0.0
    for b in b_abc_all:
        assert b > 0.0
