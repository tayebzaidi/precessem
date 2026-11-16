"""`GEO-36` step 7 (§9 item 427, §10 I24): the kept oblate spheroid in the F-human phantom.

`_f_human_kwargs(as_hole=True)` plus ``implant_spheroid = dict(SPHEROID,
keep_cells=True)`` (the F-small module's body).  Step 6's F1-F9 with the body
for the rod: (H1) CAD(403) = kept mass = body mass at STRIP_BAND and body mass
= V at CAD_ANALYTIC_BAND; (H2) CAD(3) = V_cyl - body at STRIP_BAND; (H3)
CAD(3)+CAD(403) = V_cyl at STRIP_BAND; (H4) OCC wall = A at
CAD_ANALYTIC_BAND, meshed dS(402)/A in TERMINAL_AREA_BAND; (H5) the coil
cavity did not swallow the body; (H6) census; (H7) centroid; (H8) rim
fidelity; (H9) cells < 600 000 before any facet reading.  No solve.
"""

from __future__ import annotations

import math
import time

import dolfinx
import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.io.mesh import (
    BIRDCAGE_IMPLANT_CELL_TAG,
    BIRDCAGE_IMPLANT_SURFACE_TAG,
    MeshGenerator,
)
from tests.mesh.test_birdcage_conductor_hole import (
    CAD_ANALYTIC_BAND,
    EXACT,
    _owned_census,
)
from tests.mesh.test_birdcage_implant_rod import _interior_facet_check
from tests.mesh.test_birdcage_implant_rod_f_human import (
    F_HUMAN_KEPT_ROD_CELL_GUARD,
    HOLE_CAVITY_AREA_RECORD_M2,
    HOLE_CAVITY_SURFACES_RECORD,
    HOLE_PHANTOM_CAD_RECORD_M3,
)
from tests.mesh.test_birdcage_implant_spheroid import (
    SPHEROID,
    SPHEROID_CENTRE,
    SPHEROID_SEMI_MAJOR,
    _rim_census,
    _spheroid_area,
    _spheroid_volume,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
from tests.validation.test_birdcage_f_human_rung import _ring_ports
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_th14_birdcage_copper import _outer_box_tags
from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD
from tests.validation.test_wf12_f_human_hole import _census_checks
from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs

#: First record (W1b, ``-n 1``), printed again by W2b; not asserted here.
#: Measured 417 653 = HOLE_CELL_RECORD 369 387 + 48 266 (tag 403: 14 047),
#: inside the predicted + [45 k, 115 k]
#: (`20261013T202620Z_GEO-36-step7-W1b.log:6254-6255`).
F_HUMAN_KEPT_SPHEROID_CELL_RECORD = 417653

PHANTOM_CELL_TAG = 3


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _kwargs():
    kwargs = _f_human_kwargs(as_hole=True)
    kwargs["implant_spheroid"] = dict(SPHEROID, keep_cells=True)
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
        f"\n[GEO-36 step 7] F-human kept spheroid: {n_cells} cells (guard "
        f"{F_HUMAN_KEPT_ROD_CELL_GUARD}), mesh_wall_time_s="
        f"{diag['mesh_wall_time_s']:.2f}, build elapsed {elapsed:.2f} s"
    )
    assert n_cells < F_HUMAN_KEPT_ROD_CELL_GUARD  # (H9)
    out["census"] = _owned_census(msh, cells, tdim, comm)
    (out["n402"], out["n402_ext"], out["n402_int"], out["n402_pair"],
     out["area402"], out["moments"]) = _interior_facet_check(
        msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm
    )
    one = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(1.0))
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cells)
    out["implant_volume"] = float(np.real(comm.allreduce(
        dolfinx.fem.assemble_scalar(dolfinx.fem.form(one * dx(BIRDCAGE_IMPLANT_CELL_TAG))),
        op=MPI.SUM,
    )))
    out["rim"] = _rim_census(msh, cells, comm)
    _, out["outer"] = _outer_box_tags(msh, facets, comm)
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    n_edges = comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM)
    n_faces = comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM)
    record = (
        "unregistered" if F_HUMAN_KEPT_SPHEROID_CELL_RECORD is None
        else f"{n_cells / F_HUMAN_KEPT_SPHEROID_CELL_RECORD:.6f} of "
        f"{F_HUMAN_KEPT_SPHEROID_CELL_RECORD}"
    )
    _say(
        f"[GEO-36 step 7] cells vs HOLE_CELL_RECORD {HOLE_CELL_RECORD}: difference "
        f"{n_cells - HOLE_CELL_RECORD} (predicted + [45 k, 115 k]); "
        f"F_HUMAN_KEPT_SPHEROID_CELL_RECORD {record}"
        f"\n[GEO-36 step 7] census {out['census']}; tag-403 cells "
        f"{out['census'].get(BIRDCAGE_IMPLANT_CELL_TAG, 0)}"
        f"\n[GEO-36 step 7] meshed 403 volume {out['implant_volume']:.9e} / V "
        f"{_spheroid_volume():.9e} = {out['implant_volume'] / _spheroid_volume():.6f}"
        f"\n[GEO-36 step 7] E {n_edges}, F {n_faces}; degree-1 N1curl dofs {n_edges}, "
        f"degree-2 dofs {2 * n_edges + 2 * n_faces}"
        f"\n[GEO-36 step 7] outer/cavity census {out['outer']}"
        f"\n[GEO-36 step 7] rim census {out['rim']} (n_mask, n_inside, min s, n403, "
        f"n403_outside, max s403)"
        f"\n[GEO-36 step 7] implant_spheroid_cad {diag['implant_spheroid_cad']}"
        f"\n[GEO-36 step 7] conductor_cavity_cad {diag['conductor_cavity_cad']}"
        f"\n[GEO-36 step 7] CAD masses {diag['cad_mass_by_group']}"
    )
    return out


def test_h1_implant_volume(kept):
    cad = kept["diag"]["cad_mass_by_group"]
    sph = kept["diag"]["implant_spheroid_cad"]
    body = sph["body_mass_m3"]
    h1a = abs(cad["implant"] / body - 1.0)
    h1b = abs(sph["kept_mass_m3"] / body - 1.0)
    h1c = abs(body / _spheroid_volume() - 1.0)
    _say(
        f"[GEO-36 step 7] (H1) CAD(403)/body rel {h1a:.3e}, kept/body rel {h1b:.3e} "
        f"(band {STRIP_BAND:.0e}); body {body:.12e} vs V {_spheroid_volume():.12e} rel "
        f"{h1c:.3e} (band {CAD_ANALYTIC_BAND:.0e})"
    )
    assert h1a < STRIP_BAND
    assert h1b < STRIP_BAND
    assert h1c < CAD_ANALYTIC_BAND
    assert kept["diag"]["implant_rod_cad"] is None


def test_h2_h3_partition(kept):
    cad = kept["diag"]["cad_mass_by_group"]
    body = kept["diag"]["implant_spheroid_cad"]["body_mass_m3"]
    v_cyl = _v_cyl()
    h2 = abs(cad["phantom"] / (v_cyl - body) - 1.0)
    h3 = abs((cad["phantom"] + cad["implant"]) / v_cyl - 1.0)
    _say(
        f"[GEO-36 step 7] (H2) CAD(3) {cad['phantom']:.12e} vs V_cyl - body: rel "
        f"{h2:.3e}; (H3) CAD(3)+CAD(403) vs V_cyl {v_cyl:.12e}: rel {h3:.3e} (band "
        f"{STRIP_BAND:.0e}); rod-free record {HOLE_PHANTOM_CAD_RECORD_M3:.9e} - CAD(3) "
        f"= {HOLE_PHANTOM_CAD_RECORD_M3 - cad['phantom']:.9e} vs V {_spheroid_volume():.9e}"
    )
    assert h2 < STRIP_BAND
    assert h3 < STRIP_BAND


def test_h4_wall_area(kept):
    occ = kept["diag"]["implant_spheroid_cad"]["wall_area_m2"]
    h4 = abs(occ / _spheroid_area() - 1.0)
    ratio = kept["area402"] / _spheroid_area()
    low, high = TERMINAL_AREA_BAND
    _say(
        f"[GEO-36 step 7] (H4) OCC wall {occ:.12e} vs A {_spheroid_area():.12e}: rel "
        f"{h4:.3e} (band {CAD_ANALYTIC_BAND:.0e}); meshed dS(402)/A {ratio:.6f} (band "
        f"{TERMINAL_AREA_BAND})"
    )
    assert h4 < CAD_ANALYTIC_BAND
    assert low < ratio < high


def test_h5_cavity_did_not_swallow_the_body(kept):
    cav = kept["diag"]["conductor_cavity_cad"]
    rel = abs(cav["area_m2"] / HOLE_CAVITY_AREA_RECORD_M2 - 1.0)
    _say(
        f"[GEO-36 step 7] (H5) cavity surfaces {cav['surfaces']} vs "
        f"{HOLE_CAVITY_SURFACES_RECORD}; area {cav['area_m2']:.10e} vs "
        f"{HOLE_CAVITY_AREA_RECORD_M2:.9e}: rel {rel:.3e} (band {STRIP_BAND:.0e})"
        f"\n[GEO-36 step 7] (C4s, predicted 1.181e-3) a swallowed wall misses by "
        f"{_spheroid_area() / HOLE_CAVITY_AREA_RECORD_M2:.4e}"
    )
    assert cav["surfaces"] == HOLE_CAVITY_SURFACES_RECORD
    assert rel <= STRIP_BAND


def test_h6_census(kept):
    fails = _census_checks(kept["census"], kept["n_cells"], _ring_ports())
    _say(
        f"[GEO-36 step 7] (H6) census failures {fails}; 402 owned {kept['n402']} "
        f"({kept['n402_ext']} exterior, {kept['n402_int']} interior, "
        f"{kept['n402_pair']} between tag 3 and 403)"
    )
    assert fails == []
    assert kept["census"].get(BIRDCAGE_IMPLANT_CELL_TAG, 0) > 0
    assert kept["n402"] > 0
    assert kept["n402_ext"] == 0
    assert kept["n402_int"] == kept["n402"]
    assert kept["n402_pair"] == kept["n402"]


def test_h7_wall_centroid(kept):
    c = kept["diag"]["implant_spheroid_cad"]["wall_centroid_m"]
    dev = max(abs(c[i] - SPHEROID_CENTRE[i]) for i in range(3))
    mesh_c = tuple(m / kept["area402"] for m in kept["moments"])
    _say(
        f"[GEO-36 step 7] (H7) OCC wall centroid {c} vs {SPHEROID_CENTRE}: max dev "
        f"{dev:.3e} (bound EXACT*a {EXACT * SPHEROID_SEMI_MAJOR:.1e}); meshed {mesh_c}"
    )
    assert dev <= EXACT * SPHEROID_SEMI_MAJOR


def test_h8_rim_fidelity(kept):
    n_mask, n_in, s_min, n403, n_out, s_max = kept["rim"]
    _say(
        f"[GEO-36 step 7] (H8) {n_mask} tag-3 midpoints near the rim, {n_in} inside "
        f"(min s {s_min:.6f}); tag-403 midpoints {n403}, outside {n_out} (max s "
        f"{s_max:.9f})"
    )
    assert n_mask > 0
    assert n_in == 0
    assert n403 > 0
    assert n_out == 0
