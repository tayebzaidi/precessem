"""`IMP-3` step 11 (§9 item 626, §10 I39): the rod's SAR heating the F-small gel phantom.

Fixture: `IMP-3`'s W1 record meshes (without rod 116 085 / with rod
``ROD_CELL_RECORD_W1`` = 167 371 cells) at 64 MHz, both arms through
``test_imp3_rod_birdcage._arm`` (``run_n_port_sparameter_sweep(…,
keep_fields=True)``).  The source is the **ccw superposition** (the coil's
circular drive; review reading 1 of item 626).  Inside ``_arm``'s
``field_hook``, before the fields are freed, q = G·(DG0 cell average of
½σ|E_ccw|² on tag 3) = G·SAR·ρ is moved to the tag-3 submesh, with **one**
G = 1 W / P_phantom(without rod) for both arms, so the two arms compare the
same coil drive.

Thermal side: `MAT-15` step 1's scheme imported by name — ``_cn`` (θ = ½,
N = 40, two BE half-steps), ``_be_one`` at ``DT_1``, ``_e_prime`` (the energy
helper), t_R / t* — with the gel constants of that module.  ``_run_heated`` is
bound to `MAT-15`'s 10 MHz sweep, so the pieces below it are called:
``_tag_submesh``, ``_cell_average_dg0_tag``, ``_cn``, ``_be_one``,
``_e_prime``, ``_max_imag``.

The rod wall is adiabatic.  `MAT-10`'s ``_robin_tags`` would tag every exterior
facet of the submesh, the rod's hole wall (parent tag 402) included; here Robin
(``ROBIN_TAG``) goes only on exterior facets farther than ``WALL_SELECT_M``
from the rod's capsule axis segment (the phantom's outer cylinder is ≥ 14 mm
from the rod wall, the wall itself 1 mm from the axis).  The excluded facet
count is asserted equal to the parent's tag-402 count (owned, reduced).

Asserted (every band imported): (B) ∫_sub q = G·P_arm ≤ ``B_BAND``;
(E) every run's balance ≤ ``E_BAND``; (E′) stored + surface = G·P_arm ≤
``E_BAND``; (S2) first-step slope ≤ ``S2_BAND`` — all on both arms;
(C4) the without-rod arm's 10 g θ(t*) ball means at the four 90°-rotated
centres (r = ``CENTRE_RADIUS_M``, φ = 0 on the rod's azimuth) ≤
``C4_COVARIANCE_BAND``; the census of both arms within ``CELL_COUNT_BAND``.
Printed control (predicted to miss, never asserted): the with-rod (C4).
Printed record: the tip 1 g ball θ(t*) ratio rod / without beside `IMP-3`'s
1 g SAR ratio Q on the same mesh, the whole-phantom ΔT_max ratio, P_phantom
per arm.

Scope: degree-1 EM, P1 thermal, one mesh, a 1 W normalisation, no perfusion.
Not: a published comparand ((7c), W3), compliance, F-human / F-clinical rods.
"""

from __future__ import annotations

import gc
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_SURFACE_TAG
from fem_em_solver.post.sar import averaging_ball_radius, mean_sar

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_conductor_hole import SOLID_CELL_RECORD
from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, _build_rod
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import _build
from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND
from tests.validation.test_birdcage_sar_mass_averaged import (
    CENTRE_RADIUS_M,
    PHANTOM_RHO_KG_PER_M3,
)
from tests.validation.test_imp3_rod_birdcage import (
    ROD_CELL_RECORD_W1,
    TIPS,
    _arm,
)
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    ROBIN_TAG,
    _allsum,
    _to_sub,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    DT_1,
    E_BAND,
    H_GEL,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_INF,
    T_R,
    T_STAR,
    _be_one,
    _cell_average_dg0_tag,
    _cn,
    _e_prime,
    _max_imag,
    _tag_submesh,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

KNOB = "FEM_EM_IMP3_I39"
SOURCE = "ccw"
# Capsule axis segment of `GEO-36`'s rod: 30 mm between cap centres, r = 1 mm.
ROD_HALF_SEGMENT_M = 0.015
WALL_SELECT_M = 0.005  # wall at 1 mm; nearest outer phantom surface at >= 15 mm


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _dist_to_rod_axis(x):
    """Distance of points x (3, n) from the rod's capsule axis segment."""
    zc = np.clip(x[2] - ROD_CENTRE[2], -ROD_HALF_SEGMENT_M, ROD_HALF_SEGMENT_M)
    return np.sqrt((x[0] - ROD_CENTRE[0]) ** 2 + (x[1] - ROD_CENTRE[1]) ** 2
                   + (x[2] - ROD_CENTRE[2] - zc) ** 2)


def _outer_robin_tags(sub):
    """Robin on exterior facets away from the rod wall; returns tags + reduced counts."""
    tdim = sub.topology.dim
    fdim = tdim - 1
    sub.topology.create_connectivity(fdim, tdim)
    ext = np.sort(np.asarray(dmesh.exterior_facet_indices(sub.topology), dtype=np.int32))
    mids = dmesh.compute_midpoints(sub, fdim, ext) if ext.size else np.zeros((0, 3))
    keep = _dist_to_rod_axis(mids.T) > WALL_SELECT_M if ext.size else np.zeros(0, bool)
    outer = ext[keep].astype(np.int32)
    wall = ext[~keep]
    n_owned = sub.topology.index_map(fdim).size_local
    comm = sub.comm
    counts = {
        "exterior": comm.allreduce(int(np.count_nonzero(ext < n_owned)), op=MPI.SUM),
        "outer": comm.allreduce(int(np.count_nonzero(outer < n_owned)), op=MPI.SUM),
        "wall": comm.allreduce(int(np.count_nonzero(wall < n_owned)), op=MPI.SUM),
    }
    tags = dmesh.meshtags(sub, fdim, outer, np.full(outer.size, ROBIN_TAG, dtype=np.int32))
    return tags, counts


def _parent_wall_count(msh, facet_tags):
    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    f = np.asarray(facet_tags.find(BIRDCAGE_IMPLANT_SURFACE_TAG))
    return msh.comm.allreduce(int(np.count_nonzero(f < n_owned)), op=MPI.SUM)


def _heat_arm(label, rod, msh, cell_tags, e_ccw, sigma_field, G, parent_wall):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    P = float(np.real(mean_sar(e_ccw, sigma=sigma_field, rho=PHANTOM_RHO_KG_PER_M3,
                               cell_tags=cell_tags, comm=comm,
                               subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]))
    if G is None:
        G = 1.0 / P
    dep = G * P
    sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)
    ft, counts = _outer_robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    dsall = ufl.Measure("ds", domain=sub)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    area_all = float(np.real(_allsum(comm, one * dsall)))
    Wp, q_par = _cell_average_dg0_tag(
        msh, cell_tags, PHANTOM_CELL_TAG,
        G * 0.5 * sigma_field * ufl.inner(e_ccw, e_ccw))
    q = _to_sub(Wp, q_par, sub, parent_ids)
    int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
    err_B = abs(int_q - dep) / abs(dep)
    _say(f"\nI39 [{label}] source {SOURCE}: P_phantom = {P:.12e} W; G = {G:.12e} 1/W; "
         f"G*P = {dep:.12e} W; sub cells {n_sub}; V_sub {v_sub:.12e} m^3")
    _say(f"I39 [{label}] Robin facets: exterior {counts['exterior']} = outer {counts['outer']} "
         f"+ rod-wall {counts['wall']} (parent tag-402 owned {parent_wall}); Robin area "
         f"{area:.9e} m^2, all-exterior area {area_all:.9e} m^2 (wall area "
         f"{area_all - area:.6e}, adiabatic)")
    _say(f"I39 [{label}] (B) int_sub q = {int_q:.15e} vs G*P {dep:.15e}: rel {err_B:.3e} "
         f"(bound {B_BAND:.0e})")

    t1 = time.perf_counter()
    res = _cn(sub, ft, q)
    t_cn = time.perf_counter() - t1
    ep, stor, conv = _e_prime(sub, ds, res, dep)
    res1 = _be_one(sub, ft, q)
    thbar1 = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
    ratio = RHO_C * thbar1 / (DT_1 * dep / v_sub)
    err_S2 = abs(ratio - 1.0)
    imag = max(_max_imag(res.temperature), _max_imag(res1.temperature))

    V1 = res.temperature.function_space
    th = fem.Function(V1)
    th.x.array[:] = res.temperature.x.array - T_INF
    th.x.scatter_forward()
    n_own = V1.dofmap.index_map.size_local
    loc = float(np.max(np.real(th.x.array[:n_own]))) if n_own else -np.inf
    th_max = comm.allreduce(loc, op=MPI.MAX)
    mean_th = float(np.real(_allsum(comm, th * ufl.dx(domain=sub)))) / v_sub

    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    rad_c = fem.Constant(sub, default_scalar_type(0.0))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), ufl.real(rad_c) ** 2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    f_vol = fem.form(inside * dxq)
    f_th = fem.form(th * inside * dxq)

    def _ball(centre, radius):
        centre_c.value[:] = centre
        rad_c.value = radius
        vol = float(np.real(comm.allreduce(fem.assemble_scalar(f_vol), op=MPI.SUM)))
        val = float(np.real(comm.allreduce(fem.assemble_scalar(f_th), op=MPI.SUM)))
        return val / vol, vol

    a1 = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    tip_th = {tip: _ball(np.asarray(c, float), a1) for tip, c in TIPS.items()}
    c4 = {}
    for k in range(4):
        phi = 0.5 * np.pi * k
        c4[k] = _ball(np.array([CENTRE_RADIUS_M * np.cos(phi), CENTRE_RADIUS_M * np.sin(phi), 0.0]),
                      a10)[0]
    pairs = {k: abs(c4[(k + 1) % 4] - c4[k]) / abs(c4[k]) for k in range(4)}

    _say(f"I39 [{label}] CN N={N_MAIN} to t* = {T_STAR:.3f} s (t_R {T_R:.3f} s): mean theta "
         f"{mean_th:.9e} K, max theta {th_max:.9e} K; (E) balance CN {res.max_balance_residual:.3e} "
         f"BE {res1.max_balance_residual:.3e}; max|Im theta| {imag:.3e}; CN {t_cn:.1f} s")
    _say(f"I39 [{label}] (E') stored {stor:.12e} W + surface {conv:.12e} W vs G*P {dep:.12e} W: "
         f"rel {ep:.3e} (bound {E_BAND:.0e})")
    _say(f"I39 [{label}] (S2) rho c mean(theta^1)/(dt_1 G P/V) = {ratio:.12f}, |ratio-1| "
         f"{err_S2:.4e} (bound {S2_BAND:.0e})")
    for tip, (v, vol) in tip_th.items():
        _say(f"I39 [{label}] 1 g ball (a {a1 * 1e3:.4f} mm) at {tip} tip: theta(t*) = {v:.9e} K, "
             f"ball volume in phantom {vol:.6e} m^3")
    for k in range(4):
        _say(f"I39 [{label}] (C4) 10 g theta(t*) at phi={90 * k} deg: {c4[k]:.9e} K; pair "
             f"{k}->{(k + 1) % 4} {pairs[k] * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%)")
    _say(f"I39 [{label}] thermal elapsed {time.perf_counter() - t0:.1f} s")
    return {
        "G": G, "P": P, "B": err_B, "Ep": ep, "S2": err_S2,
        "E": max(res.max_balance_residual, res1.max_balance_residual),
        "counts": counts, "parent_wall": parent_wall, "tip_th": {t: v for t, (v, _) in tip_th.items()},
        "th_max": th_max, "mean_th": mean_th, "pairs": pairs, "imag": imag,
    }


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    heat = {}

    def hook_for(parent_wall):
        def _hook(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
            G = heat["off"]["G"] if rod else None
            heat["rod" if rod else "off"] = _heat_arm(
                label, rod, msh, cell_tags, e_of[SOURCE], sigma_field, G, parent_wall)
        return _hook

    built_off = _build(True)
    off_wall = _parent_wall_count(built_off[0], built_off[2])
    off = _arm(built_off, "I39 without rod", rod=False, field_hook=hook_for(off_wall))
    del built_off
    gc.collect()
    built_rod = _build_rod()
    rod_wall = _parent_wall_count(built_rod[0], built_rod[2])
    rod = _arm(built_rod, "I39 with rod", rod=True, field_hook=hook_for(rod_wall))

    ho, hr = heat["off"], heat["rod"]
    q_sar = {t: rod["sar"][(t, SOURCE)] / off["sar"][(t, SOURCE)] for t in TIPS}
    q_th = {t: hr["tip_th"][t] / ho["tip_th"][t] for t in TIPS}
    for t in TIPS:
        _say(f"I39 record {t} tip, {SOURCE}: 1 g theta(t*) rod/without = {q_th[t]:.6f}; "
             f"IMP-3 1 g SAR ratio Q = {q_sar[t]:.6f}; ball masses without "
             f"{off['mass'][t] * 1e3:.6f} g / with {rod['mass'][t] * 1e3:.6f} g")
    _say(f"I39 record: whole-phantom max theta rod/without = {hr['th_max'] / ho['th_max']:.6f} "
         f"({hr['th_max']:.6e} / {ho['th_max']:.6e} K); mean theta rod/without "
         f"{hr['mean_th'] / ho['mean_th']:.6f}; P_phantom without {ho['P']:.9e} W, with "
         f"{hr['P']:.9e} W (dP/P {(hr['P'] - ho['P']) / ho['P']:+.4%})")
    _say(f"I39 census: without {off['cells']} (record {SOLID_CELL_RECORD}), with {rod['cells']} "
         f"(record {ROD_CELL_RECORD_W1}); arms {off['t_arm']:.1f} / {rod['t_arm']:.1f} s; "
         f"total {time.perf_counter() - t0:.1f} s at -n {comm.size}")
    return {"off": off, "rod": rod, "heat": heat}


@pytest.fixture(scope="module")
def i39():
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB}=1 not set (heavy window, ≈ 6-10 min at -n 2)")
    return _run()


# --- coding tests (W0) -------------------------------------------------------

def test_knob_refuses_without_value():
    assert os.environ.get(KNOB) in (None, "1"), f"{KNOB} must be unset or 1"


def test_wall_selector_geometry():
    pts = np.array([
        [0.016, 0.0, 0.0],        # rod side wall, 1 mm from the axis
        [0.015, 0.0, 0.016],      # top tip
        [0.015, 0.0, -0.016],     # bottom tip
        [0.03, 0.0, 0.0],         # nearest outer cylinder point
        [0.0, 0.0, 0.04],         # top end cap
        [0.015, 0.0, 0.04],       # end cap above the rod
    ]).T
    d = _dist_to_rod_axis(pts)
    np.testing.assert_allclose(d[:4], [0.001, 0.001, 0.001, 0.015], atol=1e-15)
    assert np.all(d[:3] < WALL_SELECT_M) and np.all(d[3:] > WALL_SELECT_M), d


def test_outer_robin_tags_on_a_box_keep_every_facet():
    box = dmesh.create_box(MPI.COMM_WORLD, [[-0.03, -0.03, -0.04], [0.03, 0.03, 0.04]],
                           [4, 4, 4], dmesh.CellType.tetrahedron)
    _tags, c = _outer_robin_tags(box)
    assert c["wall"] == 0 and c["outer"] == c["exterior"] == 6 * 2 * 16, c


def test_imported_scheme_is_mat15s():
    assert T_STAR == pytest.approx(0.2 * T_R, rel=1e-15) and N_MAIN == 40
    assert H_GEL == 10.0 and T_INF == 20.0


# --- W1 ------------------------------------------------------------------------

@complex_only
def test_census_both_arms(i39):
    assert abs(i39["off"]["cells"] / SOLID_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    assert abs(i39["rod"]["cells"] / ROD_CELL_RECORD_W1 - 1.0) <= CELL_COUNT_BAND


@complex_only
def test_rod_wall_is_excluded_from_robin(i39):
    off, rod = i39["heat"]["off"], i39["heat"]["rod"]
    assert off["counts"]["wall"] == 0 and off["parent_wall"] == 0, off["counts"]
    assert rod["parent_wall"] > 0
    assert rod["counts"]["wall"] == rod["parent_wall"], (rod["counts"], rod["parent_wall"])


@complex_only
@pytest.mark.parametrize("arm", ["off", "rod"])
def test_b_source_transfer(i39, arm):
    r = i39["heat"][arm]["B"]
    assert r <= B_BAND, f"(B) {arm}: {r:.3e}"


@complex_only
@pytest.mark.parametrize("arm", ["off", "rod"])
def test_e_balance(i39, arm):
    r = i39["heat"][arm]["E"]
    assert r <= E_BAND, f"(E) {arm}: {r:.3e}"


@complex_only
@pytest.mark.parametrize("arm", ["off", "rod"])
def test_e_prime_stored_plus_surface_equals_deposited(i39, arm):
    r = i39["heat"][arm]["Ep"]
    assert r <= E_BAND, f"(E') {arm}: {r:.3e}"


@complex_only
@pytest.mark.parametrize("arm", ["off", "rod"])
def test_s2_first_step_slope(i39, arm):
    r = i39["heat"][arm]["S2"]
    assert r <= S2_BAND, f"(S2) {arm}: {r:.4e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_c4_without_rod(i39, k):
    r = i39["heat"]["off"]["pairs"][k]
    assert r <= C4_COVARIANCE_BAND, f"(C4) without rod pair {k}: {r * 100:.4f}%"
