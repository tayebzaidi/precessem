"""Palace case writer and S reader (`ARCH-1` step 1a, §9 item 517, U16a).

The writer, the config dict and the reader are lifted verbatim from
``scripts/probes/ops61_step1a_palace_ans4.py`` (``export`` lines 213-302,
``config()`` 315-361, ``_read_port_s`` 400-425); ``src/`` imports neither
``scripts/`` nor ``tests/``, so these are copies pinned by equality anchors in
``tests/package/test_arch1_palace_backend.py``.  Nothing here imports or links
Palace: the run itself is a separate host-side window in ``fem-em-palace``.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

PALACE_NARROW_BASE = 20  # N_i -> 20 + i
PALACE_REMAINDER_BASE = 30  # F_i \ N_i -> 30 + i
PALACE_EXTERIOR_ATTR = 1
CASE_FILES = ("case.json", "mesh.msh")


class PalaceRunRequired(RuntimeError):
    """The case is written; Palace must run it before the read stage."""

    def __init__(self, argv, cwd):
        self.argv = list(argv)
        self.cwd = Path(cwd)
        super().__init__(f"palace run required: (cd {self.cwd} && {' '.join(self.argv)})")


@dataclass
class PalaceStudyResult:
    s_matrix: np.ndarray
    frequency_hz: float
    order: int
    cells: int
    case_dir: Path
    sha256: dict = field(default_factory=dict)
    palace_total_s: float | None = None
    backend: str = "palace"


def _tri_areas(x, tris):
    a, b, c = x[tris[:, 0]], x[tris[:, 1]], x[tris[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def write_gmsh22(path, msh, cell_tags, tri_blocks):
    """Gmsh 2.2 ASCII, the probe's operation order (lines 245-273)."""
    tdim = msh.topology.dim
    x = np.asarray(msh.geometry.x, dtype=float)
    cells = np.asarray(msh.geometry.dofmap, dtype=np.int64).reshape(-1, 4)
    ncl = msh.topology.index_map(tdim).size_local
    cells = cells[:ncl]
    ctag = np.zeros(ncl, dtype=np.int64)
    ctag[np.asarray(cell_tags.indices)] = np.asarray(cell_tags.values)
    if not np.all(ctag > 0):
        raise ValueError("an untagged tet")
    p = x[cells]
    vol6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
    neg = vol6 < 0
    cells[neg, 2], cells[neg, 3] = cells[neg, 3].copy(), cells[neg, 2].copy()
    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(path, "w") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fh.write(f"$Nodes\n{x.shape[0]}\n")
        for k, (a, b, c) in enumerate(x, start=1):
            fh.write(f"{k} {a:.17g} {b:.17g} {c:.17g}\n")
        fh.write("$EndNodes\n")
        fh.write(f"$Elements\n{ntri + ncl}\n")
        eid = 1
        for nodes, attr in tri_blocks:
            for t in nodes:
                fh.write(f"{eid} 2 2 {attr} {attr} {t[0]+1} {t[1]+1} {t[2]+1}\n")
                eid += 1
        for t, a in zip(cells, ctag):
            fh.write(f"{eid} 4 2 {a} {a} {t[0]+1} {t[1]+1} {t[2]+1} {t[3]+1}\n")
            eid += 1
        fh.write("$EndElements\n")
    return ctag


def palace_config(
    *, mesh_name, output, materials, ports, order, freq_ghz, r_ohm, conductivity=None
):
    """The probe's ``config()`` (lines 336-361) with Order / Freq as arguments.

    ``conductivity`` (``ARCH-1`` step 1f, §9 item 631): a list of Palace
    ``Boundaries.Conductivity`` entries, appended after ``LumpedPort`` as the
    probe's copper hole configs carry it
    (``scripts/probes/ops61_step1a_palace_ans4.py:1533-1536``). PEC stays on
    attribute 1 only.
    """
    cfg = {
        "Problem": {"Type": "Driven", "Verbose": 2, "Output": output},
        "Model": {"Mesh": mesh_name, "L0": 1.0},
        "Domains": {"Materials": materials},
        "Boundaries": {
            "PEC": {"Attributes": [PALACE_EXTERIOR_ATTR]},
            "LumpedPort": [
                {
                    "Index": i,
                    "R": r_ohm,
                    "Excitation": i,
                    "Attributes": [PALACE_NARROW_BASE + i],
                    "Direction": "+Z",
                }
                for i in ports
            ],
        },
        "Solver": {
            "Order": order,
            "Device": "CPU",
            "Driven": {"Samples": [{"Type": "Point", "Freq": [freq_ghz]}]},
            "Linear": {"Type": "MUMPS", "KSPType": "GMRES", "Tol": 1.0e-10, "MaxIts": 100},
        },
    }
    if conductivity is not None:
        cfg["Boundaries"]["Conductivity"] = list(conductivity)
    return cfg


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_case(
    geometry,
    built,
    case_dir,
    *,
    materials_phantom,
    order,
    freq_ghz,
    r_ohm,
    comm,
    hole_wall_sigma=None,
):
    """Write ``case.json`` / ``mesh.msh`` / ``case.sha256`` (serial only).

    Existing files are compared byte for byte; a difference raises
    ``FileExistsError`` naming the file. Never overwrites.
    """
    import dolfinx

    from ..ports._sheet_geometry import _sheet_areas

    if comm.size != 1:
        raise ValueError("the Palace case writer is serial: run it at mpiexec -n 1")
    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    msh, cell_tags, tags_n = built["mesh"], built["cell_tags"], built["facet_tags"]
    tdim = msh.topology.dim
    fdim = tdim - 1
    ports = list(range(1, geometry.leg_count + 1))
    tags_full, _areas = _sheet_areas(
        msh,
        cell_tags,
        ports,
        comm,
        geometry.sheet_iface,
        geometry.port_lower,
        geometry.port_upper,
    )

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    conductivity = None
    if hole_wall_sigma is None:
        tri_blocks = [(facet_nodes(ext), PALACE_EXTERIOR_ATTR)]
    else:
        # ARCH-1 step 1f (§9 item 631): the copper hole. Box -> attribute 1,
        # cavity wall -> its own attribute (the wall tag, as the probe's
        # ``_hole_wall_attr`` imports it, ops61_step1a_palace_ans4.py:785-788),
        # split from the exterior as ``export_hole_pec`` does (:1106, :1194-1198).
        from ..io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG

        wall_attr = int(BIRDCAGE_CONDUCTOR_SURFACE_TAG)
        # The wall tag lives on the generator's facet tags (study.py
        # ``_conductor_setup`` reads it there too), not on the sheet tags.
        wall = np.asarray(built["mesh_facet_tags"].find(wall_attr), dtype=np.int32)
        if not wall.size or not np.all(np.isin(wall, ext)):
            raise RuntimeError(f"hole wall (tag {wall_attr}): empty or not exterior")
        box = np.setdiff1d(ext, wall).astype(np.int32)
        tri_blocks = [(facet_nodes(box), PALACE_EXTERIOR_ATTR), (facet_nodes(wall), wall_attr)]
        # The Conductivity block, as ops61_step1a_palace_ans4.py:1533-1536 writes it.
        conductivity = [
            {"Attributes": [wall_attr], "Conductivity": float(hole_wall_sigma), "Permeability": 1.0}
        ]
    for i in ports:
        tag = geometry.sheet_iface + i
        full = np.asarray(tags_full.find(tag), dtype=np.int32)
        narrow = np.asarray(tags_n.find(tag), dtype=np.int32)
        if not np.all(np.isin(narrow, full)):
            raise RuntimeError(f"N_{i} not a subset of F_{i}")
        rest = np.setdiff1d(full, narrow)
        tri_blocks.append((facet_nodes(narrow), PALACE_NARROW_BASE + i))
        tri_blocks.append((facet_nodes(rest), PALACE_REMAINDER_BASE + i))

    ncl = msh.topology.index_map(tdim).size_local
    ctag = np.zeros(ncl, dtype=np.int64)
    ctag[np.asarray(cell_tags.indices)] = np.asarray(cell_tags.values)
    tags = sorted(int(t) for t in np.unique(ctag))
    air = [t for t in tags if t not in (geometry.conductor_tag, geometry.phantom_tag)]
    eps_r, sigma = materials_phantom
    phantom_air = [
        {
            "Attributes": [geometry.phantom_tag],
            "Permeability": 1.0,
            "Permittivity": eps_r,
            "Conductivity": sigma,
        },
        {"Attributes": air, "Permeability": 1.0, "Permittivity": 1.0},
    ]
    if hole_wall_sigma is None:
        materials = [
            {
                "Attributes": [geometry.conductor_tag],
                "Permeability": 1.0,
                "Permittivity": 1.0,
                "Conductivity": geometry.conductor_sigma_s_per_m,
            },
            *phantom_air,
        ]
    else:
        # The hole has no conductor volume: phantom + air, the probe's order
        # (ops61_step1a_palace_ans4.py, export_hole_pec step 5).
        materials = phantom_air
    cfg = palace_config(
        mesh_name="mesh.msh",
        output="postpro",
        materials=materials,
        ports=ports,
        order=order,
        freq_ghz=freq_ghz,
        r_ohm=r_ohm,
        conductivity=conductivity,
    )
    staging = {"case.json": json.dumps(cfg, indent=2) + "\n"}
    tmp_mesh = case_dir / "mesh.msh.new"
    write_gmsh22(tmp_mesh, msh, cell_tags, tri_blocks)
    for name in CASE_FILES:
        target = case_dir / name
        new_bytes = staging[name].encode() if name in staging else tmp_mesh.read_bytes()
        if target.exists():
            if target.read_bytes() != new_bytes:
                tmp_mesh.unlink()
                raise FileExistsError(f"{target}: differs from the Study's case")
        else:
            target.write_bytes(new_bytes)
    tmp_mesh.unlink()
    sha = {name: _sha256(case_dir / name) for name in CASE_FILES}
    manifest = "".join(f"{sha[n]}  {n}\n" for n in CASE_FILES)
    mpath = case_dir / "case.sha256"
    if mpath.exists() and mpath.read_text() != manifest:
        raise FileExistsError(f"{mpath}: manifest differs")
    mpath.write_text(manifest)
    return sha


def solve_palace(study, comm, *, case_dir, palace_ranks=8):
    """Write (or verify) the case; raise ``PalaceRunRequired`` or read S back."""
    from ..study import FourPortBirdcageGeometry, LeontovichWall, LumpedPortDrive

    geo = study.geometry
    if type(geo) is not FourPortBirdcageGeometry:
        raise TypeError(
            f"backend='palace' supports FourPortBirdcageGeometry only, " f"got {type(geo).__name__}"
        )
    if geo.implant_rod is not None:  # `WF-15` step 10: the case carries no rod-wall boundary
        raise ValueError("backend='palace': implant_rod is not supported")
    # ARCH-1 step 1f (§9 item 631): as_hole=True is accepted together with a
    # LeontovichWall conductor only (Palace Boundaries.Conductivity on the wall).
    hole_wall_sigma = None
    if geo.as_hole:
        if type(study.conductor) is not LeontovichWall:
            raise ValueError(
                "backend='palace': as_hole=True is not supported "
                "without a LeontovichWall conductor"
            )
        hole_wall_sigma = float(study.conductor.sigma_s_per_m)
    elif study.conductor is not None:
        raise ValueError("backend='palace': a conductor is not supported")
    if type(study.excitation) is not LumpedPortDrive:
        raise TypeError(
            f"backend='palace' needs a LumpedPortDrive, got " f"{type(study.excitation).__name__}"
        )
    drv = study.excitation
    if drv.z0_ohm != drv.port_impedance_ohm:
        raise ValueError(
            "backend='palace': z0_ohm must equal port_impedance_ohm "
            "(Palace's S is referenced to R)"
        )
    if case_dir is None:
        raise ValueError("backend='palace' needs case_dir")
    if comm.size != 1:
        raise ValueError("the Palace case writer is serial: run it at mpiexec -n 1")
    if study.degree not in (1, 2):
        raise ValueError(f"backend='palace': degree {study.degree} not in (1, 2)")
    freq_ghz = study.frequency_hz / 1e9
    built = geo.build(comm)
    sha = write_case(
        geo,
        built,
        case_dir,
        materials_phantom=(study.material.epsilon_r, study.material.sigma),
        order=study.degree,
        freq_ghz=freq_ghz,
        r_ohm=drv.port_impedance_ohm,
        comm=comm,
        hole_wall_sigma=hole_wall_sigma,
    )
    case_dir = Path(case_dir)
    csv_path = case_dir / "postpro" / "port-S.csv"
    if not csv_path.exists():
        raise PalaceRunRequired(["palace", "-np", str(int(palace_ranks)), "case.json"], case_dir)
    s, f = read_port_s(csv_path)
    if f != freq_ghz:
        raise ValueError(f"port-S.csv frequency {f} GHz != case Freq {freq_ghz} GHz")
    total = None
    pj = case_dir / "postpro" / "palace.json"
    if pj.exists():
        try:
            total = float(json.loads(pj.read_text())["ElapsedTime"]["Durations"]["Total"])
        except (KeyError, TypeError, ValueError):
            total = None
    return PalaceStudyResult(
        s_matrix=s,
        frequency_hz=study.frequency_hz,
        order=study.degree,
        cells=int(built["cells"]),
        case_dir=case_dir,
        sha256=sha,
        palace_total_s=total,
    )


def read_port_s(csv_path):
    """4x4 complex S and the frequency (GHz) from Palace's ``port-S.csv``."""
    with open(csv_path) as fh:
        rows = list(csv.reader(fh))
    header = [h.strip() for h in rows[0]]
    data = [float(v) for v in rows[1]]
    mag, arg = {}, {}
    freq = None
    for h, v in zip(header, data):
        m = re.search(r"S\[(\d+)\]\[(\d+)\]", h)
        if not m:
            if h.startswith("f"):
                freq = v
            continue
        key = (int(m.group(1)), int(m.group(2)))
        if "dB" in h:
            mag[key] = 10.0 ** (v / 20.0)
        elif "arg" in h:
            arg[key] = np.radians(v)
        else:
            raise ValueError(f"unparsed column {h!r}")
    n = max(max(k) for k in mag)
    s = np.full((n, n), np.nan + 0j)
    for key, a in mag.items():
        s[key[0] - 1, key[1] - 1] = a * np.exp(1j * arg[key])
    return s, freq
