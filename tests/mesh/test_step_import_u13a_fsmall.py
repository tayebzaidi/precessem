"""`WF-16` U13a′ (§9 item 1000; item 993 re-posed): the named-sheet route on the F-small four-port.

U12′ (`GEO-34` step 1d, ``tests/mesh/test_step_import_cq_port_sheet.py``)
proved on a unit box that a CadQuery-named free face reaches
``io/step_import.import_step`` by name (``highest_dim_only: false`` +
``surface_groups``), is embedded by ``occ.fragment`` and tagged with its area
exact.  This module moves that route to the F-small four-port.

Item 993 (W1a ``20261115T210648Z_WF-16-U13a-W1-n1.log:102``) read "not
embedded in any volume": its solids split each port box into two halves whose
shared face carries the strip.  Ruling (a): the fixture's solids are the
generator's **whole-box** export (``_build(False, step_export_path=...)``,
W0a′ ``20261116T030547Z_WF-16-U13a-W0a-positions-whole.log``), so each f = 0.5
strip is a face embedded inside one port box.  CadQuery assembles them with
four explicit polygon faces ``PortSheet1`` … ``PortSheet4`` at the
generator's exact-rectangle positions (W0b′
``20261116T030656Z_WF-16-U13a-W0b-author-whole.log``).  The STEPs are build
products under the gitignored ``/workspace/logs/u13a_prime/``, never
committed; the bytes read are pinned by env ``U13A_STEP_SHA256`` (W0b′'s
``[author] wrote …_whole.step … sha256=`` line).  Config
``tests/mesh/data/u13a_fsmall_four_port.json`` (step 1a's conductor / air /
phantom verbatim, one whole-box group ``port_P<i>`` = 100 + i per leg; surface
tags 230 + i).

Gates, every count and integral reduced before asserting, κ off, no solve:
(P0) ghosts; (C0) each ``port_P<i>`` claims exactly one volume; (C1) owned
strip area / (0.5·dx·dz) − 1 and / ``cad_area`` − 1 within ``STRIP_BAND``, dx,
dz from the primitive build's diagnostics; (C2′) owned facet count ==
``n_triangles``, every strip facet has exactly two cells and both carry cell
tag 100 + i (embedded in one volume); (C3) ``dS(230 + i)`` within
``STRIP_BAND·A``; (C4) ``step_name`` / ``route`` record the name-to-fragment
lineage and the global facet tag set is {1, 231, 232, 233, 234}; (C5) sheet i's
bbox-centre azimuth is 90°·(i − 1) within ``AZIMUTH_MATCH_DEG``.
Negative controls (asserted): (n-drop), (n-swap), (n-hDO), (n-absent),
(n-split) — the halves solids with step 1a's volume groups raise "is not
embedded in any volume".
(Q) printed only: min dihedral, max aspect, cells on the import beside the
primitive ``leg_sheet_strip_fraction = 0.5`` mesh (117 206-cell record).
"""

import hashlib
import os
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

import dolfinx
from dolfinx import fem

from fem_em_solver.io.step_import import import_step, load_step_config
from tests.mesh.probe_wf6_sheet_facet_census import _tet_quality
from tests.mesh.test_birdcage_leg_offset import _build
from tests.mesh.test_birdcage_port_tags import LEG_COUNT
from tests.validation.test_ops61_exact_strip_ports import (
    ANS4_EXACT_STRIP_CELL_RECORD,
    STRIP_BAND,
)
from tests.validation.test_port_birdcage_lumped_column import _sheet_bbox_centre
# `GEO-34` 1b's (iv) pre-gate band, imported from where 1b imports it
# (`test_port_birdcage_four_port` does not re-export it).
from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG

HERE = os.path.dirname(__file__)
STEP = "/workspace/logs/u13a_prime/u13a_fsmall_four_port_whole.step"
STEP_HALVES = "/workspace/logs/u13a_prime/u13a_fsmall_four_port_halves.step"
CONFIG = os.path.join(HERE, "data", "u13a_fsmall_four_port.json")
CONFIG_1A = os.path.join(HERE, "data", "geo34_fsmall_birdcage_step.json")
CELL_BASE = 100

TAG_BASE = 230
TAG_OUTER = 1
PORTS = list(range(1, LEG_COUNT + 1))

COMM = MPI.COMM_WORLD


def _p0(*args):
    if COMM.rank == 0:
        print(*args, flush=True)


def _tri_areas(mesh, facets):
    if len(facets) == 0:
        return np.zeros(0)
    g = dolfinx.mesh.entities_to_geometry(mesh, 2, np.asarray(facets, dtype=np.int32))
    x = mesh.geometry.x
    a, b, c = x[g[:, 0]], x[g[:, 1]], x[g[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def _owned(mesh, facet_tags, tag):
    n_owned = mesh.topology.index_map(2).size_local
    local = facet_tags.find(tag)
    return local[local < n_owned]


def _quality(mesh):
    """(cells, min dihedral deg, max aspect) over owned tets, reduced."""
    n_local = mesh.topology.index_map(3).size_local
    g = dolfinx.mesh.entities_to_geometry(mesh, 3, np.arange(n_local, dtype=np.int32))
    rr, dih = _tet_quality(mesh.geometry.x[g])
    lo = COMM.allreduce(float(dih.min()) if n_local else np.inf, op=MPI.MIN)
    hi = COMM.allreduce(float((1.0 / rr).max()) if n_local else 0.0, op=MPI.MAX)
    n_bad = COMM.allreduce(int(np.count_nonzero(dih < 1.0)), op=MPI.SUM)
    return int(mesh.topology.index_map(3).size_global), lo, hi, n_bad


def _azimuths(mesh, facet_tags):
    out = {}
    for i in PORTS:
        c = _sheet_bbox_centre(mesh, facet_tags, TAG_BASE + i, COMM)
        out[i] = float(np.degrees(np.arctan2(c[1], c[0])) % 360.0)
    return out


def _azimuth_error_deg(measured, target):
    return abs((measured - target + 180.0) % 360.0 - 180.0)


def _cyclic_misses(az):
    return {i: a for i, a in az.items()
            if _azimuth_error_deg(a, 90.0 * (i - 1)) >= AZIMUTH_MATCH_DEG}


@pytest.fixture(scope="module")
def primitive():
    t0 = time.perf_counter()
    msh, _cells, _facets, diag, _el = _build(
        np.zeros(LEG_COUNT), leg_sheet_strip_fraction=0.5
    )
    t_build = time.perf_counter() - t0
    dx, dy, dz = (float(v) for v in diag["port_box_size_m"])
    t0 = time.perf_counter()
    q = _quality(msh)
    t_q = time.perf_counter() - t0
    return {"dx": dx, "dy": dy, "dz": dz, "quality": q, "t_build": t_build, "t_q": t_q}


@pytest.fixture(scope="module")
def built():
    t0 = time.perf_counter()
    mesh, cell_tags, facet_tags, diag = import_step(STEP, CONFIG, COMM)
    t_import = time.perf_counter() - t0
    mesh.topology.create_connectivity(2, 3)
    _p0(f"[U13a] n={COMM.size} import + tag {t_import:.2f} s (gmsh import "
        f"{diag['import_wall_time_s']:.2f} s, mesh {diag['mesh_wall_time_s']:.2f} s)")
    return mesh, cell_tags, facet_tags, diag


def test_custody_and_config():
    pinned = os.environ["U13A_STEP_SHA256"]  # missing => KeyError: an error, not a skip
    with open(STEP, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    _p0(f"[U13a] n={COMM.size} {STEP} bytes={os.path.getsize(STEP)} sha256={digest} "
        f"U13A_STEP_SHA256={pinned}")
    assert digest == pinned
    cfg = load_step_config(CONFIG)
    assert cfg["highest_dim_only"] is False
    assert [g["step_name"] for g in cfg["surface_groups"]] == [f"PortSheet{i}" for i in PORTS]


def test_p0_ghosts(built):
    ghosts = COMM.allreduce(built[0].topology.index_map(3).num_ghosts, op=MPI.SUM)
    _p0(f"[U13a] (P0) n={COMM.size} ghost cells={ghosts}")
    assert (ghosts == 0) if COMM.size == 1 else (ghosts > 0)


def test_c0_one_volume_per_port_group(built, primitive):
    _, _, _, diag = built
    vol = primitive["dx"] * primitive["dy"] * primitive["dz"]
    for i in PORTS:
        name = f"port_P{i}"
        n_vol = diag["volumes_by_group"][name]
        _p0(f"[U13a] (C0) n={COMM.size} {name} tag={diag['tag_by_group'][name]} "
            f"volumes={n_vol} cad_mass={diag['cad_mass_by_group'][name]!r} "
            f"(dx*dy*dz={vol!r})")
        assert diag["tag_by_group"][name] == CELL_BASE + i
        assert n_vol == 1


def test_c1_strip_area(built, primitive):
    mesh, _, facet_tags, diag = built
    half = 0.5 * primitive["dx"] * primitive["dz"]
    for i in PORTS:
        rec = diag["surface_groups"][f"port_{i}"]
        a = COMM.allreduce(
            float(_tri_areas(mesh, _owned(mesh, facet_tags, TAG_BASE + i)).sum()), op=MPI.SUM
        )
        r1, r2 = a / half - 1.0, a / rec["cad_area"] - 1.0
        _p0(f"[U13a] (C1) n={COMM.size} P{i} A={a!r} A/(0.5 dx dz)-1={r1:+.3e} "
            f"A/cad_area-1={r2:+.3e} cad_area={rec['cad_area']!r} "
            f"n_triangles={rec['n_triangles']}")
        assert abs(r1) <= STRIP_BAND, (i, r1)
        assert abs(r2) <= STRIP_BAND, (i, r2)


def test_c2p_embedded_in_one_volume(built):
    mesh, cell_tags, facet_tags, diag = built
    f2c = mesh.topology.connectivity(2, 3)
    # Cell tags carried to ghosts through a DG0 scatter, so both cells of an
    # owned facet read their tag on every rank.
    V0 = fem.functionspace(mesh, ("DG", 0))
    tagf = fem.Function(V0)
    tagf.x.array[:] = -1.0
    for c, t in zip(cell_tags.indices, cell_tags.values):
        tagf.x.array[V0.dofmap.cell_dofs(int(c))[0]] = float(t)
    tagf.x.scatter_forward()
    vals = np.real(tagf.x.array)
    for i in PORTS:
        owned = _owned(mesh, facet_tags, TAG_BASE + i)
        count = COMM.allreduce(len(owned), op=MPI.SUM)
        not_two = 0
        wrong_tag = 0
        for f in owned:
            cells = f2c.links(f)
            if len(cells) != 2:
                not_two += 1
                continue
            tags = [int(round(vals[V0.dofmap.cell_dofs(int(c))[0]])) for c in cells]
            wrong_tag += int(any(t != CELL_BASE + i for t in tags))
        not_two = COMM.allreduce(not_two, op=MPI.SUM)
        wrong_tag = COMM.allreduce(wrong_tag, op=MPI.SUM)
        n_tri = diag["surface_groups"][f"port_{i}"]["n_triangles"]
        _p0(f"[U13a] (C2') n={COMM.size} P{i} owned {TAG_BASE + i}={count} "
            f"n_triangles={n_tri} facets not two-sided={not_two} "
            f"facets with a cell not tagged {CELL_BASE + i}={wrong_tag} "
            f"(embedded in one volume: 2 cells/facet, both {CELL_BASE + i})")
        assert count == n_tri
        assert not_two == 0
        assert wrong_tag == 0


def test_c3_dS(built, primitive):
    mesh, _, facet_tags, _ = built
    half = 0.5 * primitive["dx"] * primitive["dz"]
    dS = ufl.Measure("dS", domain=mesh, subdomain_data=facet_tags,
                     metadata={"quadrature_degree": 2})
    one = fem.Constant(mesh, dolfinx.default_scalar_type(1.0))
    for i in PORTS:
        v = COMM.allreduce(fem.assemble_scalar(fem.form(one * dS(TAG_BASE + i))), op=MPI.SUM)
        _p0(f"[U13a] (C3) n={COMM.size} P{i} dS({TAG_BASE + i})={v!r} rel={abs(v) / half - 1.0:+.3e}")
        assert abs(v - half) <= STRIP_BAND * half


def test_c4_lineage_and_tag_set(built):
    mesh, _, facet_tags, diag = built
    for i in PORTS:
        rec = diag["surface_groups"][f"port_{i}"]
        _p0(f"[U13a] (C4) n={COMM.size} P{i} tag={rec['tag']} step_name={rec['step_name']} "
            f"route={rec['route']}")
        assert rec["tag"] == TAG_BASE + i and rec["role"] == "port"
        assert rec["step_name"] == f"PortSheet{i}"
        assert "fragment child" in rec["route"] and "(embedded)" in rec["route"]
    n_owned = mesh.topology.index_map(2).size_local
    local = {int(t) for t in np.unique(facet_tags.values[facet_tags.indices < n_owned])}
    tags = set().union(*COMM.allgather(local))
    _p0(f"[U13a] (C4) n={COMM.size} global facet tag set={sorted(tags)}")
    assert tags == {TAG_OUTER} | {TAG_BASE + i for i in PORTS}


def test_c5_cyclic_azimuths(built):
    az = _azimuths(built[0], built[2])
    for i in PORTS:
        _p0(f"[U13a] (C5) n={COMM.size} P{i} azimuth {az[i]:.12f} deg err "
            f"{_azimuth_error_deg(az[i], 90.0 * (i - 1)):.3e} (band {AZIMUTH_MATCH_DEG:.0e})")
    assert not _cyclic_misses(az), az


def test_q_quality_printed(built, primitive):
    t0 = time.perf_counter()
    cells, dih, asp, bad = _quality(built[0])
    t_q = time.perf_counter() - t0
    pc, pd, pa, pb = primitive["quality"]
    _p0(f"[U13a] (Q) n={COMM.size} definitions: min dihedral = min over tets of the "
        f"smallest of 6 dihedral angles [deg]; max aspect = max over tets of "
        f"1/(3 r_in/R_circ) (1 = regular tet); tets with a dihedral < 1 deg counted")
    _p0(f"[U13a] (Q) n={COMM.size} import (un-optimised, no Netgen): cells={cells} "
        f"min dihedral={dih:.4f} deg max aspect={asp:.4f} tets<1deg={bad} ({t_q:.2f} s)")
    _p0(f"[U13a] (Q) n={COMM.size} primitive f=0.5 (Netgen): cells={pc} "
        f"(record {ANS4_EXACT_STRIP_CELL_RECORD}) min dihedral={pd:.4f} deg "
        f"max aspect={pa:.4f} tets<1deg={pb} (build {primitive['t_build']:.2f} s, "
        f"quality {primitive['t_q']:.2f} s)")


def test_n_drop_one_facet_is_seen(built, primitive):
    mesh, _, facet_tags, _ = built
    half = 0.5 * primitive["dx"] * primitive["dz"]
    areas = _tri_areas(mesh, _owned(mesh, facet_tags, TAG_BASE + 1))
    drop = COMM.allgather(len(areas) > 0).index(True)
    kept = areas[1:] if COMM.rank == drop else areas
    r = COMM.allreduce(float(kept.sum()), op=MPI.SUM) / half - 1.0
    _p0(f"[U13a] (n-drop) n={COMM.size} P1 drop rank {drop}: rel={r:+.3e} "
        f"= {abs(r) / STRIP_BAND:.3e} x STRIP_BAND")
    assert abs(r) > STRIP_BAND


def test_n_swap_p2_p3_fails_c5():
    cfg = load_step_config(CONFIG)
    groups = {g["name"]: g for g in cfg["surface_groups"]}
    groups["port_2"]["step_name"], groups["port_3"]["step_name"] = "PortSheet3", "PortSheet2"
    mesh, _, facet_tags, _ = import_step(STEP, cfg, COMM)
    az = _azimuths(mesh, facet_tags)
    misses = _cyclic_misses(az)
    _p0(f"[U13a] (n-swap) n={COMM.size} azimuths "
        + ", ".join(f"P{i}={v:.6f}" for i, v in sorted(az.items()))
        + f"; misses {sorted(misses)}")
    assert set(misses) == {2, 3}


def test_n_hdo_default_read_drops_free_faces():
    cfg = load_step_config(CONFIG)
    del cfg["highest_dim_only"]
    with pytest.raises(ValueError, match="no dim-2 entity named 'PortSheet1'"):
        import_step(STEP, cfg, COMM)


def test_n_absent_name_raises_everywhere():
    cfg = load_step_config(CONFIG)
    cfg["surface_groups"] = [dict(cfg["surface_groups"][0], step_name="PortSheet9")]
    with pytest.raises(ValueError, match="no dim-2 entity named 'PortSheet9'"):
        import_step(STEP, cfg, COMM)


def test_n_split_halves_interface_is_not_embedded():
    """993's interface reproduced: same four strips, port boxes split in halves."""
    cfg = load_step_config(CONFIG)
    cfg["volume_groups"] = load_step_config(CONFIG_1A)["volume_groups"]
    with open(STEP_HALVES, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    _p0(f"[U13a] (n-split) n={COMM.size} {STEP_HALVES} sha256={digest}")
    with pytest.raises(ValueError, match="is not embedded in any volume") as err:
        import_step(STEP_HALVES, cfg, COMM)
    _p0(f"[U13a] (n-split) n={COMM.size} raised: {err.value}")
