"""`GEO-36` step 2 (§9 item 171): a capped rod cut out of the F-small phantom.

``MeshGenerator.birdcage_port_domain(..., implant_rod={...})`` on the four-port
F-small fixture `tests/mesh/test_birdcage_port_sheets._build(True)` builds
(called here directly with the same imported constants, so `_build` is not
edited).  The rod: radius 1 mm, a 30 mm cylindrical section between the two
hemispherical cap centres (tip to tip 32 mm), axis z, centre (0.015, 0, 0) on
leg 1's azimuth.  r = 1 mm is the smallest rung of the step 1 probe with min
SICN >= 0.22 (`20261002T003452Z_GEO-36-probe-W1.log:49`).

Asserted (bands imported, never restated):

    (A1) CAD phantom = pi R^2 H - (pi r^2 l + 4 pi r^3 / 3) at CAD_ANALYTIC_BAND;
    (A2) sum CAD(rod) + V_rod,occ = sum CAD(flag off) at EXACT;
    (A3) OCC rod wall = 2 pi r l + 4 pi r^2 at CAD_ANALYTIC_BAND;
    (A4) meshed tag-402 area / (2 pi r l + 4 pi r^2) in TERMINAL_AREA_BAND;
    (A5) the tag census (cell set unchanged, owned census sums, 402 exterior
         with its cell in tag 3, 402 absent flag-off, flag-off = 116 085);
    (A6) OCC wall centroid on the mirror planes y = 0, z = 0 and at x = 0.015
         to EXACT * l;
    (A7) every refusal raises before any gmsh state exists.

Negative controls: (C1) the flag-off phantom misses A1's comparand by
V_rod / (pi R^2 H - V_rod) -- asserted as a miss > CAD_ANALYTIC_BAND, the
separation factor predicted and printed; (C2) predicted, printed: a fragment
that refilled the hole with air would miss A2 by V_rod / V_box; (C3) predicted,
printed: |x_bar| against EXACT * l (the x-mirror is not a symmetry).

Scope: CAD identities only.  No solve, no PEC, no SAR, no `as_hole`, no F-human.

`GEO-36` step 3 (§9 item 219, `WF-9` step 1 Part A): a third build keeps the
rod's cells (``implant_rod={..., "keep_cells": True}``; cell tag
``BIRDCAGE_IMPLANT_CELL_TAG``, 402 the interior phantom/implant interface).
Asserted: (K1) kept phantom = pi R^2 H - V_rod at CAD_ANALYTIC_BAND; (K2) CAD
implant = V_rod at CAD_ANALYTIC_BAND; (K3) sum CAD(kept) = sum CAD(flag off) at
EXACT; (K4) OCC wall = 2 pi r l + 4 pi r^2 at CAD_ANALYTIC_BAND and the meshed
402 area (a ``dS`` integral, ``("+")``) / that in TERMINAL_AREA_BAND; (K5) the
census (tag set = flag-off's + {403}, owned sum, 402 all interior and every one
between one tag-3 and one implant cell, 403 absent flag-off and in the hole
build); (K6) OCC wall centroid at (0.015, 0, 0) to EXACT * l; (K7) a non-bool
``keep_cells`` is refused before gmsh.  Controls: (C1k, asserted) the hole
build misses K3 by > EXACT; (C2k, predicted, printed) dropping 403 from the
census misses the owned sum by the implant count (ceiling 100 % of it).
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

# Bands, records, fixture and helpers imported, never restated (`ANS-1`).
from tests.mesh.test_birdcage_conductor_hole import (
    CAD_ANALYTIC_BAND,
    EXACT,
    SOLID_CELL_RECORD,
    _owned_census,
)
from tests.mesh.test_birdcage_leg_gaps import LEG_GAP_LENGTH, TERMINAL_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import (
    CELL_COUNT_BAND,
    CONDUCTOR_RESOLUTION,
)
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, _build
from tests.mesh.test_birdcage_port_tags import (
    AIR_PADDING,
    COIL_LENGTH,
    LEG_COUNT,
    LEG_SPACING,
    LEG_WIDTH,
    PHANTOM_HEIGHT,
    PHANTOM_RADIUS,
    PORT_BOX_SIZE,
    RESOLUTION,
    RING_MINOR_RADIUS,
    RING_RADIUS,
)

# The rod (§9 item 171 "Restated geometry").
ROD_RADIUS = 1.0e-3
ROD_LENGTH = 0.030  # cylindrical section between the cap centres
ROD_CENTRE = (0.015, 0.0, 0.0)
ROD_AXIS = (0.0, 0.0, 1.0)
ROD = {"radius": ROD_RADIUS, "length": ROD_LENGTH, "centre": ROD_CENTRE, "axis": ROD_AXIS}

AIR_CELL_TAG = 2
PHANTOM_CELL_TAG = 3

#: `GEO-36` step 3: the kept-rod build's global cell count, a first record
#: (printed, not asserted here; `WF-9` step 1 Part B's census reads it within
#: `CELL_COUNT_BAND`).  Registered from W1 (`-n 1`) and printed again by W2
#: (`-n 2`); predicted [180 k, 230 k] in the brief.  Measured 182 371 =
#: the hole build's 167 371 + 15 000 implant cells, the rest of the census
#: identical cell for cell (`20261004T200412Z_GEO-36-step3-W1.log:5394`).
KEPT_ROD_CELL_RECORD = 182371
#: Guard checked by `_build_kept_rod` before it returns (brief: 320 000).
KEPT_ROD_CELL_GUARD = 320_000


def _rod_volume():
    return math.pi * ROD_RADIUS**2 * ROD_LENGTH + 4.0 * math.pi * ROD_RADIUS**3 / 3.0


def _rod_area():
    return 2.0 * math.pi * ROD_RADIUS * ROD_LENGTH + 4.0 * math.pi * ROD_RADIUS**2


def _phantom_minus_rod():
    return math.pi * PHANTOM_RADIUS**2 * PHANTOM_HEIGHT - _rod_volume()


def _build_rod(**overrides):
    """The four-port F-small call `_build(True)` makes, plus ``implant_rod``."""
    comm = MPI.COMM_WORLD
    kwargs = dict(
        leg_count=LEG_COUNT,
        ring_radius=RING_RADIUS,
        leg_width=LEG_WIDTH,
        leg_spacing=LEG_SPACING,
        coil_length=COIL_LENGTH,
        ring_minor_radius=RING_MINOR_RADIUS,
        phantom_radius=PHANTOM_RADIUS,
        phantom_height=PHANTOM_HEIGHT,
        port_box_size=PORT_BOX_SIZE,
        leg_gap_length=LEG_GAP_LENGTH,
        emit_port_sheets=True,
        air_padding=AIR_PADDING,
        resolution=RESOLUTION,
        conductor_resolution=CONDUCTOR_RESOLUTION,
        phantom_resolution=None,
        comm=comm,
        return_diagnostics=True,
        implant_rod=ROD,
    )
    kwargs.update(overrides)
    started = time.perf_counter()
    out = MeshGenerator.birdcage_port_domain(**kwargs)
    return (*out, time.perf_counter() - started)


def _build_kept_rod(**overrides):
    """`GEO-36` step 3: the F-small rod build with the rod's cells kept.

    Returns ``(mesh, cell_tags, facet_tags, diagnostics, elapsed_s)``; raises
    before returning if the global cell count exceeds `KEPT_ROD_CELL_GUARD`.
    """
    out = _build_rod(implant_rod=dict(ROD, keep_cells=True), **overrides)
    msh = out[0]
    tdim = msh.topology.dim
    n_cells = msh.comm.allreduce(msh.topology.index_map(tdim).size_local, op=MPI.SUM)
    if n_cells > KEPT_ROD_CELL_GUARD:
        raise RuntimeError(
            f"_build_kept_rod: {n_cells} cells exceed the {KEPT_ROD_CELL_GUARD} guard"
        )
    return out


def _cell_values_with_ghosts(msh, cell_tags):
    """Cell tag of every local cell, ghosts included (DG0 scatter, as
    `io.mesh._interface_facet_tags` does — ghost cells need not be in
    ``cell_tags``)."""
    space = dolfinx.fem.functionspace(msh, ("DG", 0))
    marker = dolfinx.fem.Function(space)
    marker.x.array[:] = 0.0
    dofs = space.dofmap.list.reshape(-1)
    marker.x.array[dofs[cell_tags.indices]] = cell_tags.values
    marker.x.scatter_forward()
    return np.rint(np.real(marker.x.array[dofs])).astype(np.int32)


def _interior_facet_check(msh, cell_tags, facet_tags, tag, comm):
    """Global (owned facets, exterior, interior, interior between exactly one
    tag-3 and one implant cell, dS area, dS x-moments) for facet group ``tag``.

    Counted per **owned** facet; with `shared_facet` ghosting both cells of an
    owned interior facet are local.  Area and moments integrate ``dS`` on the
    ``("+")`` side, so each interior facet is counted once (`ds` reads 0 on an
    interior group — the brief's trap).
    """
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    n_owned_f = msh.topology.index_map(tdim - 1).size_local
    value = _cell_values_with_ghosts(msh, cell_tags)
    facets = np.asarray(facet_tags.indices)[np.asarray(facet_tags.values) == tag]
    owned = facets[facets < n_owned_f]
    n_ext = n_int = n_pair = 0
    for f in owned:
        cells = f2c.links(int(f))
        if len(cells) == 1:
            n_ext += 1
            continue
        n_int += 1
        if sorted(int(value[c]) for c in cells) == sorted(
            [PHANTOM_CELL_TAG, BIRDCAGE_IMPLANT_CELL_TAG]
        ):
            n_pair += 1
    x = ufl.SpatialCoordinate(msh)
    dS = ufl.Measure("dS", domain=msh, subdomain_data=facet_tags, subdomain_id=(tag,))
    one = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(1.0))

    def _total(expr):
        return float(
            np.real(comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(expr)), op=MPI.SUM))
        )

    area = _total(one * dS)
    moments = [_total(x[i]("+") * dS) for i in range(3)]
    return (
        comm.allreduce(int(owned.size), op=MPI.SUM),
        comm.allreduce(n_ext, op=MPI.SUM),
        comm.allreduce(n_int, op=MPI.SUM),
        comm.allreduce(n_pair, op=MPI.SUM),
        area,
        moments,
    )


def _facet_check(msh, cell_tags, facet_tags, tag, comm):
    """Global (facets, exterior, of which on a tag-3 cell, interior, x-moments)."""
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    n_owned_f = msh.topology.index_map(tdim - 1).size_local
    n_owned_c = msh.topology.index_map(tdim).size_local
    cell_tag_of = dict(
        zip(np.asarray(cell_tags.indices).tolist(), np.asarray(cell_tags.values).tolist())
    )
    facets = np.asarray(facet_tags.indices)[np.asarray(facet_tags.values) == tag]
    n_owned = int(np.count_nonzero(facets < n_owned_f))
    # Each exterior facet has exactly one cell, owned by exactly one rank, and
    # with `shared_facet` ghosting that rank holds every neighbour of its owned
    # cells — so counting per owned *cell* is width-invariant and a local
    # one-cell link there is truly exterior.
    n_ext_owned_cell = 0
    n_interior = 0
    n_in_phantom = 0
    for f in facets:
        cells = f2c.links(int(f))
        owned_cells = [c for c in cells if c < n_owned_c]
        if not owned_cells:
            continue
        if len(cells) == 1:
            n_ext_owned_cell += 1
            if cell_tag_of.get(int(cells[0])) == PHANTOM_CELL_TAG:
                n_in_phantom += 1
        else:
            n_interior += 1
    x = ufl.SpatialCoordinate(msh)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags, subdomain_id=(tag,))
    moments = [
        float(
            np.real(
                comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(x[i] * ds)), op=MPI.SUM)
            )
        )
        for i in range(3)
    ]
    return (
        comm.allreduce(n_owned, op=MPI.SUM),
        comm.allreduce(n_ext_owned_cell, op=MPI.SUM),
        comm.allreduce(n_in_phantom, op=MPI.SUM),
        comm.allreduce(n_interior, op=MPI.SUM),
        moments,
    )


@pytest.fixture(scope="module")
def builds():
    """Flag off first (a failed gmsh build contaminates the process, `GEO-23`)."""
    comm = MPI.COMM_WORLD
    off_msh, off_cells, off_facets, off_diag, t_off = _build(emit_port_sheets=True)
    tdim = off_msh.topology.dim
    off = {
        "n_cells": comm.allreduce(off_msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "census": _owned_census(off_msh, off_cells, tdim, comm),
        "diag": off_diag,
        "t": t_off,
        "n402": comm.allreduce(
            int(np.count_nonzero(np.asarray(off_facets.values) == BIRDCAGE_IMPLANT_SURFACE_TAG)),
            op=MPI.SUM,
        ),
    }
    rod_msh, rod_cells, rod_facets, rod_diag, t_rod = _build_rod()
    rod = {
        "n_cells": comm.allreduce(rod_msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "n_verts": comm.allreduce(rod_msh.topology.index_map(0).size_local, op=MPI.SUM),
        "census": _owned_census(rod_msh, rod_cells, tdim, comm),
        "diag": rod_diag,
        "t": t_rod,
        "area402": _exterior_facet_area(rod_msh, rod_facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm),
    }
    rod["n402"], rod["n402_ext"], rod["n402_phantom"], rod["n402_int"], moments = _facet_check(
        rod_msh, rod_cells, rod_facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm
    )
    rod["mesh_centroid"] = tuple(m / rod["area402"] for m in moments)
    if comm.rank == 0:
        print(
            f"\n[GEO-36 step 2] flag off: {off['n_cells']} cells, census "
            f"{off['census']}, mesh_wall_time_s="
            f"{off['diag']['mesh_wall_time_s']:.2f}, elapsed {off['t']:.2f} s"
            f"\n[GEO-36 step 2] rod    : {rod['n_cells']} cells {rod['n_verts']} "
            f"vertices, census {rod['census']}, tag {BIRDCAGE_IMPLANT_SURFACE_TAG} "
            f"{rod['n402']} facets, mesh_wall_time_s="
            f"{rod['diag']['mesh_wall_time_s']:.2f}, elapsed {rod['t']:.2f} s"
            f" (cell count a first record, predicted ~176 k, not asserted)"
            f"\n[GEO-36 step 2] rod CAD {rod['diag']['implant_rod_cad']}"
            f"\n[GEO-36 step 2] CAD masses off {off['diag']['cad_mass_by_group']}"
            f"\n[GEO-36 step 2] CAD masses rod {rod['diag']['cad_mass_by_group']}",
            flush=True,
        )
    # `GEO-36` step 3: the kept build last — after flag-off and the hole build
    # (`GEO-23`: a failed gmsh build contaminates its process).
    kept_msh, kept_cells, kept_facets, kept_diag, t_kept = _build_kept_rod()
    kept = {
        "n_cells": comm.allreduce(kept_msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "census": _owned_census(kept_msh, kept_cells, tdim, comm),
        "diag": kept_diag,
        "t": t_kept,
    }
    (
        kept["n402"],
        kept["n402_ext"],
        kept["n402_int"],
        kept["n402_pair"],
        kept["area402"],
        moments,
    ) = _interior_facet_check(kept_msh, kept_cells, kept_facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm)
    kept["mesh_centroid"] = tuple(m / kept["area402"] for m in moments)
    one = dolfinx.fem.Constant(kept_msh, dolfinx.default_scalar_type(1.0))
    dx = ufl.Measure("dx", domain=kept_msh, subdomain_data=kept_cells)
    kept["implant_volume"] = float(
        np.real(
            comm.allreduce(
                dolfinx.fem.assemble_scalar(dolfinx.fem.form(one * dx(BIRDCAGE_IMPLANT_CELL_TAG))),
                op=MPI.SUM,
            )
        )
    )
    if comm.rank == 0:
        record = (
            "unregistered"
            if KEPT_ROD_CELL_RECORD is None
            else f"{kept['n_cells'] / KEPT_ROD_CELL_RECORD:.6f} of {KEPT_ROD_CELL_RECORD}"
        )
        print(
            f"\n[GEO-36 step 3] kept : {kept['n_cells']} cells (first record, "
            f"predicted [180 k, 230 k], guard {KEPT_ROD_CELL_GUARD}; "
            f"KEPT_ROD_CELL_RECORD {record}), census {kept['census']}, "
            f"mesh_wall_time_s={kept['diag']['mesh_wall_time_s']:.2f}, elapsed "
            f"{kept['t']:.2f} s"
            f"\n[GEO-36 step 3] kept meshed implant volume {kept['implant_volume']:.9e}"
            f" / V_rod {_rod_volume():.9e} = {kept['implant_volume'] / _rod_volume():.6f}"
            f" (predicted [0.990, 0.999])"
            f"\n[GEO-36 step 3] kept rod CAD {kept['diag']['implant_rod_cad']}"
            f"\n[GEO-36 step 3] CAD masses kept {kept['diag']['cad_mass_by_group']}",
            flush=True,
        )
    return off, rod, kept


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def test_a1_phantom_volume_and_c1_flag_off_control(builds):
    """(A1) phantom = piR^2H - V_rod; (C1) flag off misses it (asserted)."""
    off, rod, _ = builds
    target = _phantom_minus_rod()
    a1 = abs(rod["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
    c1 = abs(off["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
    predicted = _rod_volume() / target
    _say(
        f"[GEO-36 step 2] (A1) phantom CAD {rod['diag']['cad_mass_by_group']['phantom']:.12e} "
        f"vs {target:.12e}: rel {a1:.3e} (band {CAD_ANALYTIC_BAND:.0e})"
        f"\n[GEO-36 step 2] (C1) flag-off phantom rel miss {c1:.6e} (predicted "
        f"{predicted:.6e}); separation {c1 / CAD_ANALYTIC_BAND:.4e}x the band "
        f"(predicted ceiling {predicted / CAD_ANALYTIC_BAND:.4e}x)"
    )
    assert a1 < CAD_ANALYTIC_BAND
    assert c1 > CAD_ANALYTIC_BAND


def test_a2_partition_with_the_rod_added_back(builds):
    """(A2) sum CAD(rod) + V_rod,occ = sum CAD(off); (C2) printed."""
    off, rod, _ = builds
    s_rod = sum(rod["diag"]["cad_mass_by_group"].values())
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    v_occ = rod["diag"]["implant_rod_cad"]["rod_mass_m3"]
    a2 = abs((s_rod + v_occ) / s_off - 1.0)
    c2 = v_occ / s_off
    _say(
        f"[GEO-36 step 2] (A2) {s_rod:.12e} + {v_occ:.12e} vs {s_off:.12e}: rel "
        f"{a2:.3e} (band {EXACT:.0e}); V_rod,occ vs closed form "
        f"{v_occ / _rod_volume() - 1.0:+.3e}"
        f"\n[GEO-36 step 2] (C2, predicted) an air-refilled hole would miss A2 by "
        f"{c2:.4e} = {c2 / EXACT:.4e}x EXACT"
    )
    assert a2 < EXACT


def test_a3_occ_wall_area(builds):
    _, rod, _ = builds
    area = rod["diag"]["implant_rod_cad"]["wall_area_m2"]
    a3 = abs(area / _rod_area() - 1.0)
    _say(
        f"[GEO-36 step 2] (A3) OCC wall {area:.12e} vs {_rod_area():.12e}: rel "
        f"{a3:.3e} (band {CAD_ANALYTIC_BAND:.0e})"
    )
    assert a3 < CAD_ANALYTIC_BAND


def test_a4_meshed_wall_area(builds):
    _, rod, _ = builds
    ratio = rod["area402"] / _rod_area()
    low, high = TERMINAL_AREA_BAND
    _say(
        f"[GEO-36 step 2] (A4) meshed tag {BIRDCAGE_IMPLANT_SURFACE_TAG} area "
        f"{rod['area402']:.9e} / {_rod_area():.9e} = {ratio:.6f} (band "
        f"{TERMINAL_AREA_BAND}; probe 0.99625)"
    )
    assert low < ratio < high


def test_a5_tag_census(builds):
    off, rod, _ = builds
    ports = range(1, LEG_COUNT + 1)
    expected = {1, AIR_CELL_TAG, PHANTOM_CELL_TAG}
    expected |= {PORT_LOWER + i for i in ports} | {PORT_UPPER + i for i in ports}
    ratio = off["n_cells"] / SOLID_CELL_RECORD
    _say(
        f"[GEO-36 step 2] (A5) cell tags rod {sorted(rod['census'])} off "
        f"{sorted(off['census'])}; census sums rod {sum(rod['census'].values())}/"
        f"{rod['n_cells']} off {sum(off['census'].values())}/{off['n_cells']}; "
        f"tag 402 rod {rod['n402']} facets ({rod['n402_ext']} exterior, "
        f"{rod['n402_phantom']} on a tag-3 cell, {rod['n402_int']} interior), off {off['n402']}; flag-off "
        f"{off['n_cells']}/{SOLID_CELL_RECORD} = {ratio:.6f} (band {CELL_COUNT_BAND})"
    )
    assert set(rod["census"]) == expected
    assert set(off["census"]) == expected
    assert sum(rod["census"].values()) == rod["n_cells"]
    assert sum(off["census"].values()) == off["n_cells"]
    assert rod["n402"] > 0
    assert rod["n402_int"] == 0
    assert rod["n402_ext"] == rod["n402"]
    assert rod["n402_phantom"] == rod["n402"]
    assert off["n402"] == 0
    assert abs(ratio - 1.0) < CELL_COUNT_BAND


def test_a6_mirror_symmetry_and_c3(builds):
    _, rod, _ = builds
    xb, yb, zb = rod["diag"]["implant_rod_cad"]["wall_centroid_m"]
    bound = EXACT * ROD_LENGTH
    mx, my, mz = rod["mesh_centroid"]
    _say(
        f"[GEO-36 step 2] (A6) OCC wall centroid x={xb:.12e} y={yb:+.3e} "
        f"z={zb:+.3e}; |x-0.015|={abs(xb - ROD_CENTRE[0]):.3e}, bound EXACT*l="
        f"{bound:.1e}; meshed facet centroid (printed only) x={mx:.9e} "
        f"y={my:+.3e} z={mz:+.3e}"
        f"\n[GEO-36 step 2] (C3, predicted) |x_bar|={abs(xb):.6e} = "
        f"{abs(xb) / bound:.4e}x EXACT*l: the x-mirror is not a symmetry"
    )
    assert abs(yb) <= bound
    assert abs(zb) <= bound
    assert abs(xb - ROD_CENTRE[0]) <= bound


def test_k1_k2_kept_phantom_and_implant_volumes(builds):
    """(K1) kept phantom = piR^2H - V_rod; (K2) CAD implant = V_rod."""
    _, _, kept = builds
    cad = kept["diag"]["cad_mass_by_group"]
    k1 = abs(cad["phantom"] / _phantom_minus_rod() - 1.0)
    k2 = abs(cad["implant"] / _rod_volume() - 1.0)
    _say(
        f"[GEO-36 step 3] (K1) kept phantom CAD {cad['phantom']:.12e} vs "
        f"{_phantom_minus_rod():.12e}: rel {k1:.3e} (band {CAD_ANALYTIC_BAND:.0e})"
        f"\n[GEO-36 step 3] (K2) implant CAD {cad['implant']:.12e} vs V_rod "
        f"{_rod_volume():.12e}: rel {k2:.3e} (band {CAD_ANALYTIC_BAND:.0e}); "
        f"kept_mass_m3 {kept['diag']['implant_rod_cad']['kept_mass_m3']:.12e}"
    )
    assert k1 < CAD_ANALYTIC_BAND
    assert k2 < CAD_ANALYTIC_BAND


def test_k3_nothing_removed_and_c1k_hole_control(builds):
    """(K3) sum CAD(kept) = sum CAD(off) at EXACT; (C1k) the hole misses it."""
    off, rod, kept = builds
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    s_kept = sum(kept["diag"]["cad_mass_by_group"].values())
    s_hole = sum(rod["diag"]["cad_mass_by_group"].values())
    k3 = abs(s_kept / s_off - 1.0)
    c1k = abs(s_hole / s_off - 1.0)
    _say(
        f"[GEO-36 step 3] (K3) sum CAD kept {s_kept:.12e} vs off {s_off:.12e}: rel "
        f"{k3:.3e} (band {EXACT:.0e})"
        f"\n[GEO-36 step 3] (C1k) hole build misses K3 by {c1k:.6e} (predicted "
        f"V_rod/sumCAD {_rod_volume() / s_off:.6e}; GEO-36 W2 C2 8.54e-6) = "
        f"{c1k / EXACT:.4e}x EXACT"
    )
    assert k3 < EXACT
    assert c1k > EXACT


def test_k4_kept_wall_area(builds):
    """(K4) OCC wall = closed form; meshed dS(402) area / closed form in band."""
    _, _, kept = builds
    occ_area = kept["diag"]["implant_rod_cad"]["wall_area_m2"]
    k4_occ = abs(occ_area / _rod_area() - 1.0)
    ratio = kept["area402"] / _rod_area()
    low, high = TERMINAL_AREA_BAND
    _say(
        f"[GEO-36 step 3] (K4) OCC wall {occ_area:.12e} vs {_rod_area():.12e}: rel "
        f"{k4_occ:.3e} (band {CAD_ANALYTIC_BAND:.0e}); meshed dS(402) "
        f"{kept['area402']:.9e} / {_rod_area():.9e} = {ratio:.6f} (band "
        f"{TERMINAL_AREA_BAND}; predicted 0.99625, the hole build's)"
    )
    assert k4_occ < CAD_ANALYTIC_BAND
    assert low < ratio < high


def test_k5_kept_census_and_c2k(builds):
    """(K5) the census; (C2k, predicted, printed) dropping 403 misses it."""
    off, rod, kept = builds
    n403 = kept["census"].get(BIRDCAGE_IMPLANT_CELL_TAG, 0)
    without = sum(v for t, v in kept["census"].items() if t != BIRDCAGE_IMPLANT_CELL_TAG)
    c2k_miss = kept["n_cells"] - without
    _say(
        f"[GEO-36 step 3] (K5) kept tags {sorted(kept['census'])} vs off "
        f"{sorted(off['census'])}; census sum {sum(kept['census'].values())}/"
        f"{kept['n_cells']}; implant cells {n403}; tag 402 {kept['n402']} owned "
        f"facets ({kept['n402_ext']} exterior, {kept['n402_int']} interior, "
        f"{kept['n402_pair']} between one tag-3 and one implant cell); 403 in "
        f"off {BIRDCAGE_IMPLANT_CELL_TAG in off['census']}, in hole "
        f"{BIRDCAGE_IMPLANT_CELL_TAG in rod['census']}"
        f"\n[GEO-36 step 3] (C2k, predicted) census without 403 misses by "
        f"{c2k_miss} = {c2k_miss / max(n403, 1):.4f} of the implant count "
        f"(ceiling 1.0)"
    )
    assert set(kept["census"]) == set(off["census"]) | {BIRDCAGE_IMPLANT_CELL_TAG}
    assert sum(kept["census"].values()) == kept["n_cells"]
    assert n403 > 0
    assert kept["n402"] > 0
    assert kept["n402_ext"] == 0
    assert kept["n402_int"] == kept["n402"]
    assert kept["n402_pair"] == kept["n402"]
    assert BIRDCAGE_IMPLANT_CELL_TAG not in off["census"]
    assert BIRDCAGE_IMPLANT_CELL_TAG not in rod["census"]


def test_k6_kept_wall_centroid(builds):
    """(K6) OCC wall centroid at (0.015, 0, 0) to EXACT * l; meshed printed."""
    _, _, kept = builds
    xb, yb, zb = kept["diag"]["implant_rod_cad"]["wall_centroid_m"]
    mx, my, mz = kept["mesh_centroid"]
    bound = EXACT * ROD_LENGTH
    _say(
        f"[GEO-36 step 3] (K6) OCC wall centroid x={xb:.12e} y={yb:+.3e} "
        f"z={zb:+.3e}, bound EXACT*l={bound:.1e}; meshed dS centroid (printed) "
        f"x={mx:.9e} y={my:+.3e} z={mz:+.3e}"
    )
    assert abs(xb - ROD_CENTRE[0]) <= bound
    assert abs(yb) <= bound
    assert abs(zb) <= bound


@pytest.mark.parametrize(
    "overrides",
    [
        {"as_hole": True},
        {"symmetry_quadrant": True},
        {"implant_rod": dict(ROD, radius=0.0)},
        {"implant_rod": dict(ROD, radius=-1e-3)},
        {"implant_rod": dict(ROD, length=0.0)},
        {"implant_rod": dict(ROD, length=-0.01)},
        {"implant_rod": dict(ROD, centre=(0.029, 0.0, 0.0))},
        {"implant_rod": dict(ROD, length=0.079)},
        # (K7) `GEO-36` step 3: a non-bool `keep_cells` is refused before gmsh.
        {"implant_rod": dict(ROD, keep_cells="yes")},
        {"implant_rod": dict(ROD, keep_cells=1)},
    ],
    ids=[
        "as_hole",
        "symmetry_quadrant",
        "r0",
        "rneg",
        "l0",
        "lneg",
        "radial_out",
        "axial_out",
        "k7_keep_str",
        "k7_keep_int",
    ],
)
def test_a7_refusals_raise_without_meshing(overrides):
    started = time.perf_counter()
    with pytest.raises(ValueError, match="implant_rod"):
        _build_rod(**overrides)
    assert not gmsh.isInitialized()
    _say(
        f"[GEO-36 step 2] (A7) {sorted(overrides)} refused in "
        f"{time.perf_counter() - started:.3f} s"
    )
