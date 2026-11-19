"""`GEO-36` step 5 / I18a (§9 item 316): the coated rod in the phantom, no solve.

The F-small four-port birdcage with step 3's kept rod (``keep_cells=True``,
tag 403) inside a **shaft-only dielectric sleeve** of outer radius
b = 1.5 mm (``implant_rod["coat_radius"]``, new cell tag 404): the cylinder of
radius b between the two cap-centre planes, minus the rod.  The hemispherical
caps stay bare (an insulated lead with exposed tips; review ruling (1),
2026-10-08 03:00).  The sleeve geometry is `IMP-5`'s a = 1 mm, b = 1.5 mm.

Builds flag-off -> kept -> coated in one process, in that order (`GEO-23`:
a failed gmsh build contaminates its process).  Bands and records are
imported by name, never restated; ``STRIP_BAND`` comes from
`tests/validation/test_ops61_exact_strip_ports.py` (W0 checks it collects).

Asserted: (P1)-(P9) and the negative control (C1c).  Printed only: cells per
tag, meshed volume / CAD, the three meshed interface areas against
``TERMINAL_AREA_BAND``, edge lengths in 404, E / F / dofs (the (D) identity
d1 = E, d2 = 2E + 2F asserted exactly), (C2c), (C3c).
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
    BIRDCAGE_IMPLANT_COAT_CELL_TAG,
    BIRDCAGE_IMPLANT_SURFACE_TAG,
    _interface_facet_tags,
)

# Bands, records, fixture and helpers imported, never restated (`ANS-1`).
from tests.mesh.test_birdcage_conductor_hole import EXACT, _owned_census
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_GUARD,
    KEPT_ROD_CELL_RECORD,
    PHANTOM_CELL_TAG,
    ROD,
    ROD_CENTRE,
    ROD_LENGTH,
    ROD_RADIUS,
    _build_kept_rod,
    _build_rod,
    _interior_facet_check,
    _rod_area,
    _rod_volume,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

#: `IMP-5`'s sleeve outer radius b (`IMP-5.md:26-29`); a = `ROD_RADIUS`.
COAT_RADIUS = 1.5e-3
COATED_ROD = dict(ROD, keep_cells=True, coat_radius=COAT_RADIUS)

#: Facet tags used only inside this module for the meshed tag-pair areas
#: (never emitted by the generator).
_IFACE_ROD_COAT = 9403
_IFACE_COAT_GEL = 9404
_IFACE_ROD_GEL = 9405


def _coat_volume():
    """pi (b^2 - a^2) L: the shaft-only sleeve."""
    return math.pi * (COAT_RADIUS**2 - ROD_RADIUS**2) * ROD_LENGTH


def _coat_area():
    """2 pi (a + b) L + 2 pi (b^2 - a^2): inner + outer wall + two end annuli."""
    return 2.0 * math.pi * (ROD_RADIUS + COAT_RADIUS) * ROD_LENGTH + 2.0 * math.pi * (
        COAT_RADIUS**2 - ROD_RADIUS**2
    )


def _gel_volume():
    return math.pi * PHANTOM_RADIUS**2 * PHANTOM_HEIGHT - _rod_volume() - _coat_volume()


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _total(msh, expr, comm):
    return float(
        np.real(comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(expr)), op=MPI.SUM))
    )


def _summary(msh, cells, facets, diag, elapsed, comm):
    tdim = msh.topology.dim
    return {
        "n_cells": comm.allreduce(msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "census": _owned_census(msh, cells, tdim, comm),
        "diag": diag,
        "t": elapsed,
    }


def _edge_stats(msh, cell_tags, tag, comm):
    """Global (min, mean) edge length over owned cells carrying ``tag``."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.asarray(cell_tags.indices)[np.asarray(cell_tags.values) == tag]
    idx = idx[idx < n_owned]
    x = msh.geometry.x
    dmap = msh.geometry.dofmap
    lmin = np.inf
    lsum = 0.0
    count = 0
    for c in idx:
        nodes = x[dmap[int(c)][:4]]
        for i in range(4):
            for j in range(i + 1, 4):
                length = float(np.linalg.norm(nodes[i] - nodes[j]))
                lmin = min(lmin, length)
                lsum += length
                count += 1
    gmin = comm.allreduce(lmin, op=MPI.MIN)
    gsum = comm.allreduce(lsum, op=MPI.SUM)
    gcount = comm.allreduce(count, op=MPI.SUM)
    return gmin, gsum / max(gcount, 1)


@pytest.fixture(scope="module")
def builds():
    """Flag off -> kept -> coated (`GEO-23`)."""
    comm = MPI.COMM_WORLD
    off = _summary(*_build(emit_port_sheets=True), comm)
    _say(
        f"\n[GEO-36 step 5] flag off: {off['n_cells']} cells, census {off['census']}, "
        f"mesh_wall_time_s={off['diag']['mesh_wall_time_s']:.2f}, elapsed {off['t']:.2f} s"
    )
    k_msh, k_cells, k_facets, k_diag, k_t = _build_kept_rod()
    kept = _summary(k_msh, k_cells, k_facets, k_diag, k_t, comm)
    _say(
        f"[GEO-36 step 5] kept : {kept['n_cells']} cells (KEPT_ROD_CELL_RECORD "
        f"{KEPT_ROD_CELL_RECORD}, ratio {kept['n_cells'] / KEPT_ROD_CELL_RECORD:.6f}), "
        f"census {kept['census']}, mesh_wall_time_s="
        f"{kept['diag']['mesh_wall_time_s']:.2f}, elapsed {kept['t']:.2f} s"
    )
    del k_msh, k_cells, k_facets

    msh, cells, facets, diag, t = _build_rod(implant_rod=dict(COATED_ROD))
    coated = _summary(msh, cells, facets, diag, t, comm)
    if coated["n_cells"] > KEPT_ROD_CELL_GUARD:
        raise RuntimeError(
            f"coated build: {coated['n_cells']} cells exceed the " f"{KEPT_ROD_CELL_GUARD} guard"
        )
    (
        coated["n402"],
        coated["n402_ext"],
        coated["n402_int"],
        coated["n402_pair"],
        coated["area402"],
        _,
    ) = _interior_facet_check(msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm)
    pairs = _interface_facet_tags(
        msh,
        cells,
        {
            _IFACE_ROD_COAT: (BIRDCAGE_IMPLANT_CELL_TAG, BIRDCAGE_IMPLANT_COAT_CELL_TAG),
            _IFACE_COAT_GEL: (BIRDCAGE_IMPLANT_COAT_CELL_TAG, PHANTOM_CELL_TAG),
            _IFACE_ROD_GEL: (BIRDCAGE_IMPLANT_CELL_TAG, PHANTOM_CELL_TAG),
        },
    )
    one = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(1.0))
    dS = ufl.Measure("dS", domain=msh, subdomain_data=pairs)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cells)
    coated["iface"] = {
        tag: _total(msh, one * dS(tag), comm)
        for tag in (_IFACE_ROD_COAT, _IFACE_COAT_GEL, _IFACE_ROD_GEL)
    }
    coated["vol403"] = _total(msh, one * dx(BIRDCAGE_IMPLANT_CELL_TAG), comm)
    coated["vol404"] = _total(msh, one * dx(BIRDCAGE_IMPLANT_COAT_CELL_TAG), comm)
    coated["edge404"] = _edge_stats(msh, cells, BIRDCAGE_IMPLANT_COAT_CELL_TAG, comm)
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    coated["E"] = msh.topology.index_map(1).size_global
    coated["F"] = msh.topology.index_map(2).size_global
    dofs = []
    for degree in (1, 2):
        space = dolfinx.fem.functionspace(msh, ("N1curl", degree))
        dofs.append(space.dofmap.index_map.size_global * space.dofmap.index_map_bs)
    coated["d1"], coated["d2"] = dofs
    cad = coated["diag"]["implant_rod_cad"]
    _say(
        f"[GEO-36 step 5] coated: {coated['n_cells']} cells (guard "
        f"{KEPT_ROD_CELL_GUARD}; predicted 185 k-215 k), census {coated['census']}, "
        f"403 cells {coated['census'].get(BIRDCAGE_IMPLANT_CELL_TAG, 0)} (predicted "
        f"15 000 +-10 %), 404 cells {coated['census'].get(BIRDCAGE_IMPLANT_COAT_CELL_TAG, 0)}, "
        f"mesh_wall_time_s={coated['diag']['mesh_wall_time_s']:.2f}, elapsed "
        f"{coated['t']:.2f} s"
        f"\n[GEO-36 step 5] coated rod CAD {cad}"
        f"\n[GEO-36 step 5] CAD masses coated {coated['diag']['cad_mass_by_group']}"
        f"\n[GEO-36 step 5] meshed/CAD 403 {coated['vol403']:.9e}/{_rod_volume():.9e} = "
        f"{coated['vol403'] / _rod_volume():.6f} (step 3 bare rod 0.986667); 404 "
        f"{coated['vol404']:.9e}/{_coat_volume():.9e} = "
        f"{coated['vol404'] / _coat_volume():.6f}"
        f"\n[GEO-36 step 5] 404 edge min {coated['edge404'][0]:.4e} m mean "
        f"{coated['edge404'][1]:.4e} m; E {coated['E']} F {coated['F']} d1 "
        f"{coated['d1']} d2 {coated['d2']}"
    )
    return off, kept, coated


def test_p1_p2_p3_partition(builds):
    """(P1) CAD(403) = V_rod, (P2) CAD(404) = V_coat, (P3) CAD(3) = gel."""
    _, _, coated = builds
    cad = coated["diag"]["cad_mass_by_group"]
    p1 = abs(cad["implant"] / _rod_volume() - 1.0)
    p2 = abs(cad["implant_coat"] / _coat_volume() - 1.0)
    p3 = abs(cad["phantom"] / _gel_volume() - 1.0)
    _say(
        f"[GEO-36 step 5] (P1) CAD(403) {cad['implant']:.12e} vs V_rod "
        f"{_rod_volume():.12e}: rel {p1:.3e}"
        f"\n[GEO-36 step 5] (P2) CAD(404) {cad['implant_coat']:.12e} vs V_coat "
        f"{_coat_volume():.12e}: rel {p2:.3e}"
        f"\n[GEO-36 step 5] (P3) CAD(3) {cad['phantom']:.12e} vs gel "
        f"{_gel_volume():.12e}: rel {p3:.3e} (band STRIP_BAND {STRIP_BAND:.0e})"
    )
    assert p1 < STRIP_BAND
    assert p2 < STRIP_BAND
    assert p3 < STRIP_BAND


def test_p4_nothing_removed_and_others_bitwise(builds):
    """(P4) sum CAD(coated) = sum CAD(off); groups outside {3,403,404} bitwise = kept."""
    off, kept, coated = builds
    c_cad = coated["diag"]["cad_mass_by_group"]
    k_cad = kept["diag"]["cad_mass_by_group"]
    s_c = sum(c_cad.values())
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    p4 = abs(s_c / s_off - 1.0)
    skip = {"phantom", "implant", "implant_coat"}
    others_c = {k: v for k, v in c_cad.items() if k not in skip}
    others_k = {k: v for k, v in k_cad.items() if k not in skip}
    diff = {k: others_c.get(k) == others_k.get(k) for k in set(others_c) | set(others_k)}
    _say(
        f"[GEO-36 step 5] (P4) sum CAD coated {s_c:.12e} vs off {s_off:.12e}: rel "
        f"{p4:.3e} (band {STRIP_BAND:.0e}); other groups bitwise equal to kept: "
        f"{sum(diff.values())}/{len(diff)}"
    )
    assert p4 < STRIP_BAND
    assert set(others_c) == set(others_k)
    assert all(diff.values()), {
        k: (others_c.get(k), others_k.get(k)) for k, ok in diff.items() if not ok
    }


def test_p5_occ_boundary_areas(builds):
    """(P5) OCC boundary area of 403 and of 404 against the closed forms."""
    _, _, coated = builds
    cad = coated["diag"]["implant_rod_cad"]
    a403 = cad["wall_area_m2"]
    a404 = cad["coat_wall_area_m2"]
    p5a = abs(a403 / _rod_area() - 1.0)
    p5b = abs(a404 / _coat_area() - 1.0)
    _say(
        f"[GEO-36 step 5] (P5) OCC area 403 {a403:.12e} vs {_rod_area():.12e}: rel "
        f"{p5a:.3e}; 404 {a404:.12e} vs {_coat_area():.12e}: rel {p5b:.3e} "
        f"(band {STRIP_BAND:.0e})"
    )
    assert p5a < STRIP_BAND
    assert p5b < STRIP_BAND


def test_p6_census(builds):
    """(P6) tag set = kept set + {404}; census sums; 402 between 3 and 403; guard."""
    off, kept, coated = builds
    n404 = coated["census"].get(BIRDCAGE_IMPLANT_COAT_CELL_TAG, 0)
    _say(
        f"[GEO-36 step 5] (P6) coated tags {sorted(coated['census'])} vs kept "
        f"{sorted(kept['census'])}; census sum {sum(coated['census'].values())}/"
        f"{coated['n_cells']}; 404 cells {n404}; tag 402 {coated['n402']} owned "
        f"facets ({coated['n402_ext']} exterior, {coated['n402_int']} interior, "
        f"{coated['n402_pair']} between one tag-3 and one 403 cell); 404 in kept "
        f"{BIRDCAGE_IMPLANT_COAT_CELL_TAG in kept['census']}, in off "
        f"{BIRDCAGE_IMPLANT_COAT_CELL_TAG in off['census']}; cells "
        f"{coated['n_cells']} < guard {KEPT_ROD_CELL_GUARD}"
    )
    assert set(coated["census"]) == set(kept["census"]) | {BIRDCAGE_IMPLANT_COAT_CELL_TAG}
    assert sum(coated["census"].values()) == coated["n_cells"]
    assert n404 > 0
    assert BIRDCAGE_IMPLANT_COAT_CELL_TAG not in kept["census"]
    assert BIRDCAGE_IMPLANT_COAT_CELL_TAG not in off["census"]
    assert coated["n402"] > 0
    assert coated["n402_ext"] == 0
    assert coated["n402_int"] == coated["n402"]
    assert coated["n402_pair"] == coated["n402"]
    assert coated["n_cells"] < KEPT_ROD_CELL_GUARD


def test_p7_coat_centroid(builds):
    """(P7) OCC coat centroid at (0.015, 0, 0) to EXACT * l (the K6 pattern)."""
    _, _, coated = builds
    xb, yb, zb = coated["diag"]["implant_rod_cad"]["coat_centroid_m"]
    bound = EXACT * ROD_LENGTH
    _say(
        f"[GEO-36 step 5] (P7) OCC coat centroid x={xb:.12e} y={yb:+.3e} "
        f"z={zb:+.3e}; |x-0.015|={abs(xb - ROD_CENTRE[0]):.3e}, bound {bound:.1e}"
    )
    assert abs(xb - ROD_CENTRE[0]) <= bound
    assert abs(yb) <= bound
    assert abs(zb) <= bound


def test_p9_kept_record(builds):
    """(P9) the kept build in this window reads KEPT_ROD_CELL_RECORD."""
    _, kept, _ = builds
    ratio = kept["n_cells"] / KEPT_ROD_CELL_RECORD
    _say(
        f"[GEO-36 step 5] (P9) kept {kept['n_cells']} / KEPT_ROD_CELL_RECORD "
        f"{KEPT_ROD_CELL_RECORD} = {ratio:.6f} (band {CELL_COUNT_BAND})"
    )
    assert abs(ratio - 1.0) < CELL_COUNT_BAND


def test_c1c_kept_build_misses_p3(builds):
    """(C1c, asserted) the uncoated kept build misses P3 by V_coat / gel."""
    _, kept, _ = builds
    c1c = abs(kept["diag"]["cad_mass_by_group"]["phantom"] / _gel_volume() - 1.0)
    predicted = _coat_volume() / _gel_volume()
    _say(
        f"[GEO-36 step 5] (C1c) kept phantom misses P3 by {c1c:.6e} (predicted "
        f"{predicted:.6e}) = {c1c / STRIP_BAND:.4e}x STRIP_BAND"
    )
    assert c1c > STRIP_BAND


def test_c2c_c3c_printed(builds):
    """(C2c, C3c, predicted, printed) a capped sleeve; a 403/404 swap."""
    _, _, coated = builds
    cad = coated["diag"]["cad_mass_by_group"]
    capped = _coat_volume() + 4.0 * math.pi * (COAT_RADIUS**3 - ROD_RADIUS**3) / 3.0
    c2c = abs(capped / _coat_volume() - 1.0)
    c3c = abs(cad["implant_coat"] / _rod_volume() - 1.0)
    _say(
        f"[GEO-36 step 5] (C2c, predicted 8.44e-2) a capped sleeve misses P2 by "
        f"{c2c:.4e}; measured P2 miss {abs(cad['implant_coat'] / _coat_volume() - 1.0):.3e}"
        f"\n[GEO-36 step 5] (C3c, predicted 0.1968) a 403/404 swap misses P1 by {c3c:.4f}"
    )


def test_interface_areas_printed(builds):
    """Printed, never gated: meshed tag-pair dS areas / closed form."""
    _, _, coated = builds
    a, b, L = ROD_RADIUS, COAT_RADIUS, ROD_LENGTH
    refs = {
        _IFACE_ROD_COAT: ("(403|404)", 2.0 * math.pi * a * L),
        _IFACE_COAT_GEL: ("(404|3)", 2.0 * math.pi * b * L + 2.0 * math.pi * (b**2 - a**2)),
        _IFACE_ROD_GEL: ("(403|3)", 4.0 * math.pi * a**2),
    }
    for tag, (name, ref) in refs.items():
        area = coated["iface"][tag]
        _say(
            f"[GEO-36 step 5] meshed {name} {area:.9e} / {ref:.9e} = {area / ref:.6f} "
            f"(predicted inside TERMINAL_AREA_BAND {TERMINAL_AREA_BAND}, printed)"
        )
    _say(
        f"[GEO-36 step 5] meshed dS(402) {coated['area402']:.9e} (= the (403|3) "
        f"interface, rebuilt by the generator)"
    )


def test_dof_identity(builds):
    """(D) d1 = E, d2 = 2E + 2F, exact."""
    _, _, coated = builds
    _say(
        f"[GEO-36 step 5] (D) E {coated['E']} F {coated['F']} d1 {coated['d1']} d2 "
        f"{coated['d2']} = 2E+2F {2 * coated['E'] + 2 * coated['F']}"
    )
    assert coated["d1"] == coated["E"]
    assert coated["d2"] == 2 * coated["E"] + 2 * coated["F"]


@pytest.mark.parametrize(
    "coat",
    [
        dict(COATED_ROD, coat_radius=True),
        dict(COATED_ROD, coat_radius="1.5e-3"),
        dict(COATED_ROD, coat_radius=float("nan")),
        dict(COATED_ROD, coat_radius=float("inf")),
        dict(COATED_ROD, coat_radius=ROD_RADIUS),
        dict(COATED_ROD, coat_radius=0.5 * ROD_RADIUS),
        dict(ROD, coat_radius=COAT_RADIUS),
        dict(ROD, keep_cells=False, coat_radius=COAT_RADIUS),
        dict(COATED_ROD, coat_radius=PHANTOM_RADIUS),
    ],
    ids=["bool", "str", "nan", "inf", "eq_a", "lt_a", "no_keep", "keep_false", "outside"],
)
def test_p8_refusals_raise_without_meshing(coat):
    """(P8) refused before gmsh, ``gmsh.isInitialized()`` False (K7 pattern)."""
    started = time.perf_counter()
    with pytest.raises(ValueError, match="implant_rod"):
        _build_rod(implant_rod=coat)
    assert not gmsh.isInitialized()
    _say(
        f"[GEO-36 step 5] (P8) {coat.get('coat_radius')!r} / keep "
        f"{coat.get('keep_cells')!r} refused in {time.perf_counter() - started:.3f} s"
    )
