"""`GEO-36` step 6 (§9 item 372, §10 I21a): the kept rod in the F-human hole phantom.

``MeshGenerator.birdcage_port_domain(**_f_human_kwargs(as_hole=True),
implant_rod=dict(ROD, keep_cells=True))`` — the H9′ copper-hole mesh with the
F-small rod (r = 1 mm, l = 30 mm, axis z, centre (0.015, 0, 0)) kept as tag 403.
CAD identities only; no solve.  Every band imported (`ANS-1`).

Asserted: (F1) CAD(403) = V_rod and kept_mass = rod_mass = V_rod; (F2) CAD(3) =
V_cyl - V_rod; (F3) CAD(3) + CAD(403) = V_cyl; (F4) OCC wall = A_rod and meshed
dS(402) / A_rod in TERMINAL_AREA_BAND; (F5) the coil cavity did not swallow the
rod wall (240 surfaces, area at its record); (F6) census; (F7) OCC wall centroid;
(F8) refusals before gmsh; (F9) cell guard.  Controls: (C1h, asserted) the
rod-free phantom record misses F2 by > STRIP_BAND; (C2h, C3h) predicted, printed.
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
from tests.mesh.test_birdcage_conductor_hole import EXACT, _owned_census
from tests.mesh.test_birdcage_implant_rod import (  # noqa: F401  (used through _interior_facet_check)
    ROD,
    ROD_CENTRE,
    ROD_LENGTH,
    ROD_RADIUS,
    _cell_values_with_ghosts,
    _interior_facet_check,
    _rod_area,
    _rod_volume,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
from tests.validation.test_birdcage_f_human_rung import _ring_ports
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_th14_birdcage_copper import _outer_box_tags
from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs
from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD
from tests.validation.test_wf12_f_human_hole import _census_checks

ROD_FH = dict(ROD, keep_cells=True)

#: First record (W1, ``-n 1``), printed again by W2; I21b asserts its census
#: against it at `CELL_COUNT_BAND`.  Printed, not asserted here.  Measured
#: 436 541 = HOLE_CELL_RECORD 369 387 + 67 154 (tag 403: 15 025), inside the
#: predicted [420 k, 455 k] (`20261011T185757Z_GEO-36-step6-W1.log:6251`).
F_HUMAN_KEPT_ROD_CELL_RECORD = 436541
F_HUMAN_KEPT_ROD_CELL_GUARD = 600_000

# Copied at the printed digits from
# `20261010T005723Z_GEO-35-s2a-W1-wf12-buildonly.log:88, :6219` (the rod-free
# F-human hole build): ten significant figures, rounding <= 2.4e-10 relative.
HOLE_CAVITY_SURFACES_RECORD = 240
HOLE_CAVITY_AREA_RECORD_M2 = 2.094760146e-01
HOLE_PHANTOM_CAD_RECORD_M3 = 2.225676428e-03

PHANTOM_CELL_TAG = 3


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _kwargs():
    kwargs = _f_human_kwargs(as_hole=True)
    kwargs["implant_rod"] = dict(ROD_FH)
    return kwargs


def _v_cyl():
    k = _f_human_kwargs(as_hole=True)
    return math.pi * k["phantom_radius"] ** 2 * k["phantom_height"]


@pytest.fixture(scope="module")
def kept():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **_kwargs()
    )
    elapsed = time.perf_counter() - t0
    tdim = msh.topology.dim
    n_cells = comm.allreduce(msh.topology.index_map(tdim).size_local, op=MPI.SUM)
    out = {"n_cells": n_cells, "diag": diag, "t": elapsed}
    _say(
        f"\n[GEO-36 step 6] F-human kept rod: {n_cells} cells (guard "
        f"{F_HUMAN_KEPT_ROD_CELL_GUARD}), mesh_wall_time_s="
        f"{diag['mesh_wall_time_s']:.2f}, build elapsed {elapsed:.2f} s"
    )
    # (F9) before any further reading.
    assert n_cells < F_HUMAN_KEPT_ROD_CELL_GUARD
    out["census"] = _owned_census(msh, cells, tdim, comm)
    out["n402"], out["n402_ext"], out["n402_int"], out["n402_pair"], out["area402"], moments = (
        _interior_facet_check(msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm)
    )
    one = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(1.0))
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cells)
    out["implant_volume"] = float(
        np.real(
            comm.allreduce(
                dolfinx.fem.assemble_scalar(dolfinx.fem.form(one * dx(BIRDCAGE_IMPLANT_CELL_TAG))),
                op=MPI.SUM,
            )
        )
    )
    _, out["outer"] = _outer_box_tags(msh, facets, comm)
    msh.topology.create_entities(1)
    n_edges = comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM)
    n_faces = comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM)
    out["E"], out["F"] = n_edges, n_faces
    record = (
        "unregistered"
        if F_HUMAN_KEPT_ROD_CELL_RECORD is None
        else f"{n_cells / F_HUMAN_KEPT_ROD_CELL_RECORD:.6f} of {F_HUMAN_KEPT_ROD_CELL_RECORD}"
    )
    n403 = out["census"].get(BIRDCAGE_IMPLANT_CELL_TAG, 0)
    _say(
        f"[GEO-36 step 6] cells vs HOLE_CELL_RECORD {HOLE_CELL_RECORD}: ratio "
        f"{n_cells / HOLE_CELL_RECORD:.6f}, difference {n_cells - HOLE_CELL_RECORD} "
        f"(predicted [420 k, 455 k]); F_HUMAN_KEPT_ROD_CELL_RECORD {record}"
        f"\n[GEO-36 step 6] census {out['census']}; tag-403 cells {n403} "
        f"(predicted 15 000 +/- 10 %); tag-3 cells "
        f"{out['census'].get(PHANTOM_CELL_TAG, 0)}"
        f"\n[GEO-36 step 6] meshed 403 volume {out['implant_volume']:.9e} / V_rod "
        f"{_rod_volume():.9e} = {out['implant_volume'] / _rod_volume():.6f} "
        f"(F-small kept 0.986667)"
        f"\n[GEO-36 step 6] E {n_edges}, F {n_faces}; degree-1 N1curl dofs {n_edges}, "
        f"degree-2 dofs {2 * n_edges + 2 * n_faces}"
        f"\n[GEO-36 step 6] outer/cavity census {out['outer']} (401 facets vs "
        f"78 066 rod-free, printed only)"
        f"\n[GEO-36 step 6] implant_rod_cad {diag['implant_rod_cad']}"
        f"\n[GEO-36 step 6] conductor_cavity_cad {diag['conductor_cavity_cad']}"
        f"\n[GEO-36 step 6] CAD masses {diag['cad_mass_by_group']}"
    )
    out["moments"] = moments
    return out


def test_f1_implant_volume(kept):
    cad = kept["diag"]["cad_mass_by_group"]
    rod_cad = kept["diag"]["implant_rod_cad"]
    v = _rod_volume()
    f1 = abs(cad["implant"] / v - 1.0)
    f1k = abs(rod_cad["kept_mass_m3"] / v - 1.0)
    f1r = abs(rod_cad["rod_mass_m3"] / v - 1.0)
    c3h = cad["phantom"] / cad["implant"]
    _say(
        f"[GEO-36 step 6] (F1) CAD(403) {cad['implant']:.12e} vs V_rod {v:.12e}: rel "
        f"{f1:.3e}; kept_mass rel {f1k:.3e}; rod_mass rel {f1r:.3e} (band "
        f"{STRIP_BAND:.0e})\n[GEO-36 step 6] (C3h, predicted ~2.3e4) swapping 3/403 "
        f"misses F1 by {c3h:.4e}"
    )
    assert f1 < STRIP_BAND
    assert f1k < STRIP_BAND
    assert f1r < STRIP_BAND


def test_f2_f3_partition_and_c1h(kept):
    cad = kept["diag"]["cad_mass_by_group"]
    v_cyl = _v_cyl()
    target = v_cyl - _rod_volume()
    f2 = abs(cad["phantom"] / target - 1.0)
    f3 = abs((cad["phantom"] + cad["implant"]) / v_cyl - 1.0)
    c1h = abs(HOLE_PHANTOM_CAD_RECORD_M3 / target - 1.0)
    _say(
        f"[GEO-36 step 6] (F2) CAD(3) {cad['phantom']:.12e} vs V_cyl - V_rod "
        f"{target:.12e}: rel {f2:.3e}; (F3) CAD(3)+CAD(403) vs V_cyl {v_cyl:.12e}: "
        f"rel {f3:.3e} (band {STRIP_BAND:.0e})"
        f"\n[GEO-36 step 6] rod-free record {HOLE_PHANTOM_CAD_RECORD_M3:.9e} - CAD(3) = "
        f"{HOLE_PHANTOM_CAD_RECORD_M3 - cad['phantom']:.9e} vs V_rod {_rod_volume():.9e}"
        f"\n[GEO-36 step 6] (C1h) rod-free record misses F2 by {c1h:.4e} (predicted "
        f"4.42e-5) = {c1h / STRIP_BAND:.4e}x STRIP_BAND"
    )
    assert f2 < STRIP_BAND
    assert f3 < STRIP_BAND
    assert c1h > STRIP_BAND


def test_f4_wall_area(kept):
    occ = kept["diag"]["implant_rod_cad"]["wall_area_m2"]
    f4 = abs(occ / _rod_area() - 1.0)
    ratio = kept["area402"] / _rod_area()
    low, high = TERMINAL_AREA_BAND
    _say(
        f"[GEO-36 step 6] (F4) OCC wall {occ:.12e} vs A_rod {_rod_area():.12e}: rel "
        f"{f4:.3e} (band {STRIP_BAND:.0e}); meshed dS(402)/A_rod {ratio:.6f} (band "
        f"{TERMINAL_AREA_BAND}; predicted 0.99625)"
    )
    assert f4 < STRIP_BAND
    assert low < ratio < high


def test_f5_cavity_did_not_swallow_the_rod(kept):
    cav = kept["diag"]["conductor_cavity_cad"]
    rel = abs(cav["area_m2"] / HOLE_CAVITY_AREA_RECORD_M2 - 1.0)
    _say(
        f"[GEO-36 step 6] (F5) cavity surfaces {cav['surfaces']} vs "
        f"{HOLE_CAVITY_SURFACES_RECORD}; area {cav['area_m2']:.10e} vs "
        f"{HOLE_CAVITY_AREA_RECORD_M2:.9e}: rel {rel:.3e} (band {STRIP_BAND:.0e})"
        f"\n[GEO-36 step 6] (C2h, predicted) a swallowed wall misses by "
        f"{_rod_area() / HOLE_CAVITY_AREA_RECORD_M2:.4e}"
    )
    assert cav["surfaces"] == HOLE_CAVITY_SURFACES_RECORD
    assert rel <= STRIP_BAND


def test_f6_census(kept):
    fails = _census_checks(kept["census"], kept["n_cells"], _ring_ports())
    _say(
        f"[GEO-36 step 6] (F6) census failures {fails}; 402 owned {kept['n402']} "
        f"({kept['n402_ext']} exterior, {kept['n402_int']} interior, "
        f"{kept['n402_pair']} between tag 3 and 403)"
    )
    assert fails == []
    assert kept["census"].get(BIRDCAGE_IMPLANT_CELL_TAG, 0) > 0
    assert kept["n402"] > 0
    assert kept["n402_ext"] == 0
    assert kept["n402_int"] == kept["n402"]
    assert kept["n402_pair"] == kept["n402"]


def test_f7_wall_centroid(kept):
    c = kept["diag"]["implant_rod_cad"]["wall_centroid_m"]
    dev = max(abs(c[i] - ROD_CENTRE[i]) for i in range(3))
    mesh_c = tuple(m / kept["area402"] for m in kept["moments"])
    _say(
        f"[GEO-36 step 6] (F7) OCC wall centroid {c} vs {ROD_CENTRE}: max dev "
        f"{dev:.3e} (bound EXACT*l {EXACT * ROD_LENGTH:.1e}); meshed {mesh_c}"
    )
    assert dev <= EXACT * ROD_LENGTH


def _rh():
    k = _f_human_kwargs(as_hole=True)
    return k["phantom_radius"], k["phantom_height"]


@pytest.mark.parametrize(
    "rod, kwargs",
    [
        (ROD, {"as_hole": True}),
        (dict(ROD, keep_cells=True, coat_radius=1.5e-3), {"as_hole": True}),
        (ROD_FH, {"symmetry_quadrant": True}),
    ],
    ids=["hole-as_hole", "coat-as_hole", "kept-quadrant"],
)
def test_f8_refusals(rod, kwargs):
    r, h = _rh()
    with pytest.raises(ValueError, match="implant_rod"):
        MeshGenerator._normalise_implant_rod(rod, phantom_radius=r, phantom_height=h, **kwargs)
    assert not gmsh.isInitialized()


def test_f8_kept_rod_accepted_with_as_hole():
    r, h = _rh()
    out = MeshGenerator._normalise_implant_rod(
        ROD_FH, phantom_radius=r, phantom_height=h, as_hole=True
    )
    assert out["keep_cells"] is True
    assert not gmsh.isInitialized()
