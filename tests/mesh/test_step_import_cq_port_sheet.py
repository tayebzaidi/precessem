"""`GEO-34` step 1d / U12' (§9 item 985): a CadQuery-authored STEP with a named
free face, mapped to ``port`` by name and tagged natively by ``import_step``.

Fixture: ``tests/mesh/data/geo34_step1d_port_sheet.step``, authored once by
``tests/mesh/cq_author_geo34_step1d.py`` under CadQuery 2.8.0 in window W0a
(``20261114T181641Z_GEO-34-U12p-W0a-author.log:39``: 21468 bytes,
sha256 pinned below).  W0b (``20261114T181729Z_GEO-34-U12p-W0b-read.log``) read
the face as ``Shapes/geo34_step1d/PortSheet/PortSheet`` only with
``highestDimOnly=False`` (``:136``; the default read r1 has six unnamed faces and
no PortSheet, ``:60–67``), the fragment child embedded (``:184``
``getEmbedded(3,1)=[(2, 7)]``) and outside the wall set (``:182``).

Real build at any width; every count and integral reduced before asserting.
Unit box with the rectangle (0.1, 0.2, 0.47, 0.8 x 0.6), constants of
``tests/mesh/probe_geo37_interior_facet.py:30–32``.
"""

import hashlib
import os

import dolfinx
import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.io.step_import import import_step, load_step_config
from fem_em_solver.ports.definitions import (
    PortDefinition,
    validate_required_port_tags_exist,
)

HERE = os.path.dirname(__file__)
STEP = os.path.join(HERE, "data", "geo34_step1d_port_sheet.step")
CONFIG = os.path.join(HERE, "data", "geo34_step1d_port_sheet.json")
SHA256 = "f729d87d88fee04de4f44524fe0edc204f273bf00123e99a6c37fce3c3fde392"

RECT_AREA = 0.8 * 0.6
H = 0.08
BOX_AREA = 6.0
STRIP_BAND = 1.0e-9  # tests/validation/test_ops61_exact_strip_ports.py:58
TAG_PORT, TAG_WALL = 201, 101

COMM = MPI.COMM_WORLD


def _tri_areas(mesh, facets):
    if len(facets) == 0:
        return np.zeros(0)
    g = dolfinx.mesh.entities_to_geometry(mesh, 2, np.asarray(facets, dtype=np.int32))
    x = mesh.geometry.x
    a, b, c = x[g[:, 0]], x[g[:, 1]], x[g[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


@pytest.fixture(scope="module")
def built():
    mesh, cell_tags, facet_tags, diag = import_step(STEP, CONFIG, COMM)
    return mesh, cell_tags, facet_tags, diag


def _owned(mesh, facet_tags, tag):
    n_owned = mesh.topology.index_map(2).size_local
    local = facet_tags.find(tag)
    return local[local < n_owned]


def test_custody_and_config():
    with open(STEP, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    print(f"[U12p] n={COMM.size} sha256={digest}")
    assert digest == SHA256
    cfg = load_step_config(CONFIG)
    assert cfg["resolution"] == H


def test_p0_ghosts(built):
    mesh = built[0]
    ghosts = COMM.allreduce(mesh.topology.index_map(3).num_ghosts, op=MPI.SUM)
    print(f"[U12p] (P0) n={COMM.size} ghost cells={ghosts}")
    if COMM.size == 1:
        assert ghosts == 0
    else:
        assert ghosts > 0


def test_p1_cad_identity(built):
    mesh, _, facet_tags, diag = built
    rec = diag["surface_groups"]["port_1"]
    assert rec["tag"] == TAG_PORT and rec["role"] == "port"
    assert rec["step_name"] == "PortSheet"
    assert abs(rec["cad_area"] - RECT_AREA) <= 1e-12
    area = COMM.allreduce(
        float(_tri_areas(mesh, _owned(mesh, facet_tags, TAG_PORT)).sum()), op=MPI.SUM
    )
    r1 = area / RECT_AREA - 1.0
    r2 = area / rec["cad_area"] - 1.0
    print(
        f"[U12p] (P1) n={COMM.size} A={area!r} A/0.48-1={r1:+.3e} "
        f"A/cad_area-1={r2:+.3e} n_triangles={rec['n_triangles']} route={rec['route']}"
    )
    assert abs(r1) <= STRIP_BAND
    assert abs(r2) <= STRIP_BAND


def test_p2_complete_and_embedded(built):
    mesh, _, facet_tags, diag = built
    owned = _owned(mesh, facet_tags, TAG_PORT)
    count = COMM.allreduce(len(owned), op=MPI.SUM)
    mesh.topology.create_connectivity(2, 3)
    f2c = mesh.topology.connectivity(2, 3)
    one_sided = sum(1 for f in owned if len(f2c.links(f)) != 2)
    one_sided = COMM.allreduce(one_sided, op=MPI.SUM)
    print(
        f"[U12p] (P2) n={COMM.size} owned 201={count} "
        f"n_triangles={diag['surface_groups']['port_1']['n_triangles']} one-sided={one_sided}"
    )
    assert count == diag["surface_groups"]["port_1"]["n_triangles"]
    assert one_sided == 0


def test_p3_dS_and_w_instrument(built):
    mesh, _, facet_tags, _ = built
    md = {"quadrature_degree": 2}
    dS = ufl.Measure("dS", domain=mesh, subdomain_data=facet_tags, metadata=md)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags, metadata=md)
    one = fem.Constant(mesh, dolfinx.default_scalar_type(1.0))
    i_port = COMM.allreduce(fem.assemble_scalar(fem.form(one * dS(TAG_PORT))), op=MPI.SUM)
    i_wall = COMM.allreduce(fem.assemble_scalar(fem.form(one * ds(TAG_WALL))), op=MPI.SUM)
    r3 = abs(i_port) / RECT_AREA - 1.0
    rw = abs(i_wall) / BOX_AREA - 1.0
    print(
        f"[U12p] (P3) n={COMM.size} dS(201)={i_port!r} rel={r3:+.3e}; "
        f"(W) ds(101)={i_wall!r} rel={rw:+.3e}"
    )
    assert abs(i_port - RECT_AREA) <= STRIP_BAND * RECT_AREA
    assert abs(i_wall - BOX_AREA) <= STRIP_BAND * BOX_AREA


def test_p4_tag_set_and_validator(built):
    mesh, _, facet_tags, _ = built
    n_owned = mesh.topology.index_map(2).size_local
    local = {int(t) for t in np.unique(facet_tags.values[facet_tags.indices < n_owned])}
    tags = set().union(*COMM.allgather(local))
    print(f"[U12p] (P4) n={COMM.size} global facet tag set={sorted(tags)}")
    assert tags == {TAG_WALL, TAG_PORT}
    validate_required_port_tags_exist(
        [PortDefinition("P1", TAG_PORT, TAG_WALL, "cq_named_sheet")], tags
    )


def test_n3_dropped_facet_is_seen(built):
    mesh, _, facet_tags, _ = built
    areas = _tri_areas(mesh, _owned(mesh, facet_tags, TAG_PORT))
    owners = COMM.allgather(len(areas) > 0)
    drop = owners.index(True)
    kept = areas[1:] if COMM.rank == drop else areas
    a = COMM.allreduce(float(kept.sum()), op=MPI.SUM)
    r = a / RECT_AREA - 1.0
    print(
        f"[U12p] (n3) n={COMM.size} drop rank {drop}: A'/0.48-1={r:+.3e} "
        f"= {abs(r) / STRIP_BAND:.3e} x STRIP_BAND (predicted ~5e6)"
    )
    assert abs(r) > STRIP_BAND


def test_n_hdo_default_read_drops_free_face():
    """W0b r1 (``20261114T181729Z_GEO-34-U12p-W0b-read.log:60–67``): default
    ``highestDimOnly`` reads six unnamed faces and no PortSheet — asserted."""
    cfg = load_step_config(CONFIG)
    del cfg["highest_dim_only"]
    with pytest.raises(ValueError, match="PortSheet"):
        import_step(STEP, cfg, COMM)


def test_n1_absent_name_raises_everywhere():
    cfg = load_step_config(CONFIG)
    cfg["surface_groups"] = [dict(cfg["surface_groups"][0], step_name="PortSheetX")]
    with pytest.raises(ValueError, match="PortSheetX"):
        import_step(STEP, cfg, COMM)


def test_n2_without_surface_groups_no_port_tag():
    cfg = load_step_config(CONFIG)
    del cfg["surface_groups"]
    mesh, _, facet_tags, diag = import_step(STEP, cfg, COMM)
    local = {int(t) for t in np.unique(facet_tags.values)}
    tags = set().union(*COMM.allgather(local))
    print(f"[U12p] (n2) n={COMM.size} tag set without surface_groups={sorted(tags)}")
    assert TAG_PORT not in tags
    assert "surface_groups" not in diag
