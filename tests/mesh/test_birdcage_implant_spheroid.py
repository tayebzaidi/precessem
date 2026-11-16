"""`GEO-36` step 7 (§9 item 427, §10 I24): an oblate PEC spheroid in the F-small phantom.

``MeshGenerator.birdcage_port_domain(..., implant_spheroid={...})`` on the
four-port F-small fixture (`tests/mesh/test_birdcage_implant_rod._build_rod`'s
kwargs with ``implant_rod=None``).  The body: a = 6 mm equatorial, c = 1.2 mm
along x-hat (aspect 5, `IMP-6`'s FEM-gated body), centre (0.015, 0, 0); built
by revolving an exact half-ellipse, so it is the exact ellipsoid of revolution
and its rim circle is a model edge.  Rim points (0.015, 0, +-a).

Asserted (bands imported, never restated):

    (S1) OCC body mass = 4 pi a^2 c / 3 at CAD_ANALYTIC_BAND;
    (S2) OCC wall area = the oblate closed form at CAD_ANALYTIC_BAND;
    (S3) hole: CAD(3) = pi R^2 H - V at CAD_ANALYTIC_BAND; sum CAD(hole) +
         body mass = sum CAD(off) at EXACT;
    (S4) kept: CAD(403) = kept mass = body mass at STRIP_BAND; CAD(3)+CAD(403)
         = pi R^2 H at STRIP_BAND; sum CAD(kept) = sum CAD(off) at EXACT;
    (S5) meshed 402 area / A in TERMINAL_AREA_BAND (ds hole, dS("+") kept);
    (S6) the census (hole / kept / flag-off = SOLID_CELL_RECORD);
    (S7) OCC wall centroid at the centre to EXACT * a;
    (S8) rim fidelity: owned tag-3 cells with midpoint within 2 rho_r of a rim
         point exist, and none of their midpoints lies inside the true
         spheroid; kept: every tag-403 midpoint inside it;
    (S9) refusals before any gmsh state;
    (S10) cells < KEPT_ROD_CELL_GUARD before any facet reading.

Negative controls: (C1s, asserted) the `addSphere` + `dilate` body misses
(S1)'s comparand by > CAD_ANALYTIC_BAND (4.374e-4 measured, `IMP-6` step 1
W0); (C2s, asserted) the flag-off phantom misses (S3)'s comparand by
V / (pi R^2 H - V); (C3s, predicted) the hole misses (S4)'s EXACT sum.
No solve, no PEC condition, no SAR.
"""

from __future__ import annotations

import math
import time

import dolfinx
import gmsh
import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.io.mesh import (
    BIRDCAGE_IMPLANT_CELL_TAG,
    BIRDCAGE_IMPLANT_SURFACE_TAG,
    MeshGenerator,
)
from tests.mesh.probe_two_torus_conductor_hole import _exterior_facet_area
from tests.mesh.test_birdcage_conductor_hole import (
    CAD_ANALYTIC_BAND,
    EXACT,
    SOLID_CELL_RECORD,
    _owned_census,
)
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_GUARD,
    PHANTOM_CELL_TAG,
    ROD,
    _build_rod,
    _cell_values_with_ghosts,  # noqa: F401  (used through _interior_facet_check)
    _interior_facet_check,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

SPHEROID_SEMI_MAJOR = 0.006
SPHEROID_SEMI_MINOR = 0.0012
SPHEROID_CENTRE = (0.015, 0.0, 0.0)
SPHEROID_AXIS = (1.0, 0.0, 0.0)
SPHEROID = {
    "semi_major": SPHEROID_SEMI_MAJOR,
    "semi_minor": SPHEROID_SEMI_MINOR,
    "centre": SPHEROID_CENTRE,
    "axis": SPHEROID_AXIS,
}
RIM_POINTS = {"top": (0.015, 0.0, +0.006), "bottom": (0.015, 0.0, -0.006)}
#: `IMP-6` step 1 W0: the dilated-sphere body's volume miss at this a, c
#: (`20261012T184107Z_IMP-6-step1-W0.log:503`).
DILATE_MISS_RECORD = 4.374e-4

#: First records (W1a, ``-n 1``), printed again by W2a; not asserted here.
#: Measured 150 140 = 116 085 + 34 055 (predicted + [35 k, 85 k]: a reading
#: just under the range, rule (e)) and 164 266 = hole + 14 126 tag-403 cells
#: (`20261013T202258Z_GEO-36-step7-W1a.log:5390-5391`).
SPHEROID_HOLE_CELL_RECORD = 150140
KEPT_SPHEROID_CELL_RECORD = 164266


def _rho_r():
    return SPHEROID_SEMI_MINOR**2 / SPHEROID_SEMI_MAJOR


def _spheroid_volume():
    return 4.0 * math.pi * SPHEROID_SEMI_MAJOR**2 * SPHEROID_SEMI_MINOR / 3.0


def _spheroid_area():
    a, c = SPHEROID_SEMI_MAJOR, SPHEROID_SEMI_MINOR
    e = math.sqrt(1.0 - c * c / (a * a))
    return 2.0 * math.pi * a * a * (1.0 + (1.0 - e * e) / e * math.atanh(e))


def _inside_true_spheroid(points):
    """s = rho_l^2/a^2 + z_l^2/c^2 in the body frame (s < 1 inside)."""
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    d = pts - np.asarray(SPHEROID_CENTRE)
    ax = np.asarray(SPHEROID_AXIS) / np.linalg.norm(SPHEROID_AXIS)
    z_l = d @ ax
    rho2 = np.maximum(np.sum(d * d, axis=1) - z_l * z_l, 0.0)
    return rho2 / SPHEROID_SEMI_MAJOR**2 + z_l * z_l / SPHEROID_SEMI_MINOR**2


def _build_spheroid(**overrides):
    kwargs = dict(implant_rod=None, implant_spheroid=dict(SPHEROID))
    kwargs.update(overrides)
    return _build_rod(**kwargs)


def _build_kept_spheroid(**overrides):
    out = _build_spheroid(implant_spheroid=dict(SPHEROID, keep_cells=True), **overrides)
    msh = out[0]
    tdim = msh.topology.dim
    n_cells = msh.comm.allreduce(msh.topology.index_map(tdim).size_local, op=MPI.SUM)
    if n_cells > KEPT_ROD_CELL_GUARD:
        raise RuntimeError(
            f"_build_kept_spheroid: {n_cells} cells exceed the {KEPT_ROD_CELL_GUARD} guard"
        )
    return out


def _rim_census(msh, cell_tags, comm):
    """(S8): owned-cell midpoint census near the rim, reduced over ranks.

    Returns (n_mask, n_mask_inside, min s in mask, n403, n403_outside, max s 403).
    """
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.asarray(cell_tags.indices)
    val = np.asarray(cell_tags.values)
    keep = idx < n_owned
    idx, val = idx[keep], val[keep]
    rim = np.asarray(list(RIM_POINTS.values()))
    radius = 2.0 * _rho_r()
    ph = idx[val == PHANTOM_CELL_TAG].astype(np.int32)
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, ph)
    near = np.zeros(mid.shape[0], dtype=bool)
    for p in rim:
        near |= np.linalg.norm(mid - p, axis=1) < radius
    s_near = _inside_true_spheroid(mid[near])
    n_mask = int(near.sum())
    n_in = int(np.count_nonzero(s_near < 1.0))
    s_min = float(s_near.min()) if s_near.size else np.inf
    imp = idx[val == BIRDCAGE_IMPLANT_CELL_TAG].astype(np.int32)
    s_imp = _inside_true_spheroid(dolfinx.mesh.compute_midpoints(msh, tdim, imp))
    n_out = int(np.count_nonzero(s_imp >= 1.0 + 1.0e-9))
    s_max = float(s_imp.max()) if s_imp.size else -np.inf
    return (
        comm.allreduce(n_mask, op=MPI.SUM),
        comm.allreduce(n_in, op=MPI.SUM),
        comm.allreduce(s_min, op=MPI.MIN),
        comm.allreduce(int(imp.size), op=MPI.SUM),
        comm.allreduce(n_out, op=MPI.SUM),
        comm.allreduce(s_max, op=MPI.MAX),
    )


def _exterior_402_check(msh, cell_tags, facet_tags, comm):
    """Hole build: (owned 402 facets, exterior with its cell in tag 3, interior)."""
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    n_owned_f = msh.topology.index_map(tdim - 1).size_local
    value = _cell_values_with_ghosts(msh, cell_tags)
    facets = np.asarray(facet_tags.indices)[
        np.asarray(facet_tags.values) == BIRDCAGE_IMPLANT_SURFACE_TAG
    ]
    owned = facets[facets < n_owned_f]
    n_ext3 = n_int = 0
    for f in owned:
        cells = f2c.links(int(f))
        if len(cells) == 1:
            if int(value[cells[0]]) == PHANTOM_CELL_TAG:
                n_ext3 += 1
        else:
            n_int += 1
    return (
        comm.allreduce(int(owned.size), op=MPI.SUM),
        comm.allreduce(n_ext3, op=MPI.SUM),
        comm.allreduce(n_int, op=MPI.SUM),
    )


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _edges_faces(msh, comm):
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    return (
        comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM),
        comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM),
    )


@pytest.fixture(scope="module")
def builds():
    """Flag off, then hole, then kept (`GEO-23`: one order, never a retry)."""
    comm = MPI.COMM_WORLD
    off_msh, off_cells, _, off_diag, t_off = _build(emit_port_sheets=True)
    tdim = off_msh.topology.dim
    off = {
        "n_cells": comm.allreduce(off_msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "census": _owned_census(off_msh, off_cells, tdim, comm),
        "diag": off_diag,
        "t": t_off,
    }
    _say(f"\n[GEO-36 step 7] flag off: {off['n_cells']} cells, census {off['census']}, "
         f"elapsed {t_off:.2f} s")

    h_msh, h_cells, h_facets, h_diag, t_h = _build_spheroid()
    hole = {
        "n_cells": comm.allreduce(h_msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "diag": h_diag,
        "t": t_h,
    }
    _say(f"[GEO-36 step 7] hole: {hole['n_cells']} cells (guard {KEPT_ROD_CELL_GUARD}), "
         f"mesh_wall_time_s={h_diag['mesh_wall_time_s']:.2f}, elapsed {t_h:.2f} s")
    assert hole["n_cells"] < KEPT_ROD_CELL_GUARD  # (S10)
    hole["census"] = _owned_census(h_msh, h_cells, tdim, comm)
    hole["area402"] = _exterior_facet_area(h_msh, h_facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm)
    hole["n402"], hole["n402_ext3"], hole["n402_int"] = _exterior_402_check(
        h_msh, h_cells, h_facets, comm
    )
    x = ufl.SpatialCoordinate(h_msh)
    ds = ufl.Measure("ds", domain=h_msh, subdomain_data=h_facets,
                     subdomain_id=(BIRDCAGE_IMPLANT_SURFACE_TAG,))
    hole["mesh_centroid"] = tuple(
        float(np.real(comm.allreduce(
            dolfinx.fem.assemble_scalar(dolfinx.fem.form(x[i] * ds)), op=MPI.SUM
        ))) / hole["area402"]
        for i in range(3)
    )
    hole["rim"] = _rim_census(h_msh, h_cells, comm)
    hole["E"], hole["F"] = _edges_faces(h_msh, comm)

    k_msh, k_cells, k_facets, k_diag, t_k = _build_kept_spheroid()
    kept = {
        "n_cells": comm.allreduce(k_msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "diag": k_diag,
        "t": t_k,
    }
    _say(f"[GEO-36 step 7] kept: {kept['n_cells']} cells, mesh_wall_time_s="
         f"{k_diag['mesh_wall_time_s']:.2f}, elapsed {t_k:.2f} s")
    assert kept["n_cells"] < KEPT_ROD_CELL_GUARD  # (S10)
    kept["census"] = _owned_census(k_msh, k_cells, tdim, comm)
    (kept["n402"], kept["n402_ext"], kept["n402_int"], kept["n402_pair"],
     kept["area402"], moments) = _interior_facet_check(
        k_msh, k_cells, k_facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm
    )
    kept["mesh_centroid"] = tuple(m / kept["area402"] for m in moments)
    one = dolfinx.fem.Constant(k_msh, dolfinx.default_scalar_type(1.0))
    dx = ufl.Measure("dx", domain=k_msh, subdomain_data=k_cells)
    kept["implant_volume"] = float(np.real(comm.allreduce(
        dolfinx.fem.assemble_scalar(dolfinx.fem.form(one * dx(BIRDCAGE_IMPLANT_CELL_TAG))),
        op=MPI.SUM,
    )))
    kept["rim"] = _rim_census(k_msh, k_cells, comm)
    kept["E"], kept["F"] = _edges_faces(k_msh, comm)

    def _rec(n, rec):
        return "unregistered" if rec is None else f"{n / rec:.6f} of {rec}"

    _say(
        f"[GEO-36 step 7] hole cells {hole['n_cells']} = {SOLID_CELL_RECORD} + "
        f"{hole['n_cells'] - SOLID_CELL_RECORD} (predicted + [35 k, 85 k]); "
        f"SPHEROID_HOLE_CELL_RECORD {_rec(hole['n_cells'], SPHEROID_HOLE_CELL_RECORD)}"
        f"\n[GEO-36 step 7] kept cells {kept['n_cells']} = hole + "
        f"{kept['n_cells'] - hole['n_cells']}; tag-403 cells "
        f"{kept['census'].get(BIRDCAGE_IMPLANT_CELL_TAG, 0)} (predicted [10 k, 30 k]); "
        f"KEPT_SPHEROID_CELL_RECORD {_rec(kept['n_cells'], KEPT_SPHEROID_CELL_RECORD)}"
        f"\n[GEO-36 step 7] hole census {hole['census']}\n[GEO-36 step 7] kept census "
        f"{kept['census']}"
        f"\n[GEO-36 step 7] meshed 403 volume {kept['implant_volume']:.9e} / V "
        f"{_spheroid_volume():.9e} = {kept['implant_volume'] / _spheroid_volume():.6f} "
        f"(predicted [0.97, 0.995])"
        f"\n[GEO-36 step 7] hole E {hole['E']} F {hole['F']} (N1curl dofs deg1 "
        f"{hole['E']}, deg2 {2 * hole['E'] + 2 * hole['F']}); kept E {kept['E']} F "
        f"{kept['F']} (deg1 {kept['E']}, deg2 {2 * kept['E'] + 2 * kept['F']})"
        f"\n[GEO-36 step 7] rim census hole {hole['rim']} kept {kept['rim']} "
        f"(n_mask, n_inside, min s, n403, n403_outside, max s403)"
        f"\n[GEO-36 step 7] hole spheroid CAD {h_diag['implant_spheroid_cad']}"
        f"\n[GEO-36 step 7] kept spheroid CAD {k_diag['implant_spheroid_cad']}"
        f"\n[GEO-36 step 7] implant_rod_cad hole {h_diag['implant_rod_cad']} kept "
        f"{k_diag['implant_rod_cad']}"
        f"\n[GEO-36 step 7] CAD masses off {off_diag['cad_mass_by_group']}"
        f"\n[GEO-36 step 7] CAD masses hole {h_diag['cad_mass_by_group']}"
        f"\n[GEO-36 step 7] CAD masses kept {k_diag['cad_mass_by_group']}"
    )
    return off, hole, kept


def test_s1_s2_body_mass_and_wall_area(builds):
    _, hole, kept = builds
    for name, b in (("hole", hole), ("kept", kept)):
        cad = b["diag"]["implant_spheroid_cad"]
        s1 = abs(cad["body_mass_m3"] / _spheroid_volume() - 1.0)
        s2 = abs(cad["wall_area_m2"] / _spheroid_area() - 1.0)
        _say(
            f"[GEO-36 step 7] (S1) {name} body mass {cad['body_mass_m3']:.12e} vs V "
            f"{_spheroid_volume():.12e}: rel {s1:.3e}; (S2) wall {cad['wall_area_m2']:.12e}"
            f" vs A {_spheroid_area():.12e}: rel {s2:.3e} (band {CAD_ANALYTIC_BAND:.0e}); "
            f"seam bbox {cad['seam_bbox_m']}"
        )
        assert s1 < CAD_ANALYTIC_BAND
        assert s2 < CAD_ANALYTIC_BAND
        assert b["diag"]["implant_rod_cad"] is None


def test_s3_hole_phantom_and_partition_and_c2s(builds):
    off, hole, _ = builds
    v_cyl = math.pi * PHANTOM_RADIUS**2 * PHANTOM_HEIGHT
    target = v_cyl - _spheroid_volume()
    s3a = abs(hole["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    s_hole = sum(hole["diag"]["cad_mass_by_group"].values())
    body = hole["diag"]["implant_spheroid_cad"]["body_mass_m3"]
    s3b = abs((s_hole + body) / s_off - 1.0)
    c2s = abs(off["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
    pred = _spheroid_volume() / target
    _say(
        f"[GEO-36 step 7] (S3) hole CAD(3) vs piR^2H - V {target:.12e}: rel {s3a:.3e} "
        f"(band {CAD_ANALYTIC_BAND:.0e}); sum CAD(hole)+body vs off: rel {s3b:.3e} "
        f"(band {EXACT:.0e})"
        f"\n[GEO-36 step 7] (C2s) flag-off phantom misses by {c2s:.6e} (predicted "
        f"{pred:.6e}) = {c2s / CAD_ANALYTIC_BAND:.4e}x the band (ceiling "
        f"{pred / CAD_ANALYTIC_BAND:.4e}x)"
    )
    assert s3a < CAD_ANALYTIC_BAND
    assert s3b < EXACT
    assert c2s > CAD_ANALYTIC_BAND


def test_s4_kept_partition_and_c3s(builds):
    off, hole, kept = builds
    cad = kept["diag"]["cad_mass_by_group"]
    sph = kept["diag"]["implant_spheroid_cad"]
    v_cyl = math.pi * PHANTOM_RADIUS**2 * PHANTOM_HEIGHT
    a = abs(cad["implant"] / sph["body_mass_m3"] - 1.0)
    b = abs(sph["kept_mass_m3"] / sph["body_mass_m3"] - 1.0)
    c = abs((cad["phantom"] + cad["implant"]) / v_cyl - 1.0)
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    s_kept = sum(cad.values())
    s_hole = sum(hole["diag"]["cad_mass_by_group"].values())
    d = abs(s_kept / s_off - 1.0)
    c3s = abs(s_hole / s_off - 1.0)
    _say(
        f"[GEO-36 step 7] (S4) CAD(403)/body rel {a:.3e}, kept_mass/body rel {b:.3e}, "
        f"CAD(3)+CAD(403) vs piR^2H rel {c:.3e} (band {STRIP_BAND:.0e}); sum kept vs "
        f"off rel {d:.3e} (band {EXACT:.0e})"
        f"\n[GEO-36 step 7] (C3s, predicted ~1.6e-5) hole misses the sum by {c3s:.6e} "
        f"(body/sumCAD {sph['body_mass_m3'] / s_off:.6e})"
    )
    assert a < STRIP_BAND
    assert b < STRIP_BAND
    assert c < STRIP_BAND
    assert d < EXACT


def test_s5_meshed_wall_area(builds):
    _, hole, kept = builds
    low, high = TERMINAL_AREA_BAND
    rh = hole["area402"] / _spheroid_area()
    rk = kept["area402"] / _spheroid_area()
    _say(
        f"[GEO-36 step 7] (S5) meshed 402 / A: hole (ds) {rh:.6f}, kept (dS+) {rk:.6f} "
        f"(band {TERMINAL_AREA_BAND}; predicted [0.985, 0.999])"
    )
    assert low < rh < high
    assert low < rk < high


def test_s6_census(builds):
    off, hole, kept = builds
    ratio = off["n_cells"] / SOLID_CELL_RECORD
    _say(
        f"[GEO-36 step 7] (S6) tags off {sorted(off['census'])} hole "
        f"{sorted(hole['census'])} kept {sorted(kept['census'])}; sums hole "
        f"{sum(hole['census'].values())}/{hole['n_cells']} kept "
        f"{sum(kept['census'].values())}/{kept['n_cells']}; hole 402 {hole['n402']} "
        f"({hole['n402_ext3']} exterior on tag 3, {hole['n402_int']} interior); kept 402 "
        f"{kept['n402']} ({kept['n402_ext']} exterior, {kept['n402_int']} interior, "
        f"{kept['n402_pair']} 3|403); flag-off {off['n_cells']}/{SOLID_CELL_RECORD} = "
        f"{ratio:.6f} (band {CELL_COUNT_BAND})"
    )
    assert set(hole["census"]) == set(off["census"])
    assert sum(hole["census"].values()) == hole["n_cells"]
    assert hole["n402"] > 0
    assert hole["n402_int"] == 0
    assert hole["n402_ext3"] == hole["n402"]
    assert BIRDCAGE_IMPLANT_CELL_TAG not in hole["census"]
    assert set(kept["census"]) == set(off["census"]) | {BIRDCAGE_IMPLANT_CELL_TAG}
    assert sum(kept["census"].values()) == kept["n_cells"]
    assert kept["n402"] > 0
    assert kept["n402_ext"] == 0
    assert kept["n402_int"] == kept["n402"]
    assert kept["n402_pair"] == kept["n402"]
    assert abs(ratio - 1.0) < CELL_COUNT_BAND


def test_s7_wall_centroid(builds):
    _, hole, kept = builds
    bound = EXACT * SPHEROID_SEMI_MAJOR
    for name, b in (("hole", hole), ("kept", kept)):
        c = b["diag"]["implant_spheroid_cad"]["wall_centroid_m"]
        _say(
            f"[GEO-36 step 7] (S7) {name} OCC wall centroid ({c[0]:.12e}, {c[1]:+.3e}, "
            f"{c[2]:+.3e}), bound EXACT*a {bound:.1e}; meshed (printed) "
            f"{tuple(f'{v:.9e}' for v in b['mesh_centroid'])}"
        )
        for i in range(3):
            assert abs(c[i] - SPHEROID_CENTRE[i]) <= bound


def test_s8_rim_fidelity(builds):
    _, hole, kept = builds
    for name, b in (("hole", hole), ("kept", kept)):
        n_mask, n_in, s_min, n403, n_out, s_max = b["rim"]
        _say(
            f"[GEO-36 step 7] (S8) {name}: {n_mask} tag-3 midpoints within 2 rho_r of a "
            f"rim point, {n_in} inside the true spheroid (min s {s_min:.6f}); tag-403 "
            f"midpoints {n403}, outside {n_out} (max s {s_max:.9f})"
        )
        assert n_mask > 0
        assert n_in == 0
    assert kept["rim"][3] > 0
    assert kept["rim"][4] == 0


@pytest.mark.parametrize(
    "overrides",
    [
        {"implant_rod": dict(ROD)},
        {"symmetry_quadrant": True},
        {"as_hole": True},
        {"implant_spheroid": dict(SPHEROID, semi_minor=SPHEROID_SEMI_MAJOR)},
        {"implant_spheroid": dict(SPHEROID, semi_minor=0.0)},
        {"implant_spheroid": dict(SPHEROID, axis=(0.0, 0.0, 0.0))},
        {"implant_spheroid": dict(SPHEROID, centre=(0.025, 0.0, 0.0))},
        {"implant_spheroid": dict(SPHEROID, keep_cells="yes")},
    ],
    ids=["rod_and_spheroid", "symmetry_quadrant", "as_hole_hole", "a_eq_c", "c0",
         "zero_axis", "poking_out", "keep_str"],
)
def test_s9_refusals_raise_without_meshing(overrides):
    started = time.perf_counter()
    with pytest.raises(ValueError, match="implant_spheroid"):
        _build_spheroid(**overrides)
    assert not gmsh.isInitialized()
    _say(f"[GEO-36 step 7] (S9) {sorted(overrides)} refused in "
         f"{time.perf_counter() - started:.3f} s")


def test_s9_kept_accepted_with_as_hole():
    out = MeshGenerator._normalise_implant_spheroid(
        dict(SPHEROID, keep_cells=True),
        phantom_radius=PHANTOM_RADIUS,
        phantom_height=PHANTOM_HEIGHT,
        as_hole=True,
    )
    assert out["keep_cells"] is True
    assert not gmsh.isInitialized()


def test_c1s_dilated_sphere_misses_s1():
    """(C1s, asserted) `addSphere` + `dilate` misses 4 pi a^2 c / 3."""
    a, c = SPHEROID_SEMI_MAJOR, SPHEROID_SEMI_MINOR
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("c1s_dilate")
        s = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.rotate([(3, s)], 0, 0, 0, 0, 0, 1, 0.25 * math.pi)
        gmsh.model.occ.dilate([(3, s)], 0, 0, 0, a, a, c)
        gmsh.model.occ.synchronize()
        mass = float(gmsh.model.occ.getMass(3, s))
    finally:
        gmsh.finalize()
    miss = abs(mass / _spheroid_volume() - 1.0)
    _say(
        f"[GEO-36 step 7] (C1s) dilated-sphere mass {mass:.12e} vs V "
        f"{_spheroid_volume():.12e}: miss {miss:.4e} (IMP-6 record "
        f"{DILATE_MISS_RECORD:.3e}) = {miss / CAD_ANALYTIC_BAND:.4e}x the band"
    )
    assert miss > CAD_ANALYTIC_BAND
