"""`GEO-34` step 1d W0 (§9 item 289): do face / volume names survive Gmsh's
own OCC STEP writer and come back through ``occ.importShapes``?

Probe script (never collected: ``python_files = test_*.py``), one rank.  The
fixture is `GEO-37` item 222's unit box with one embedded rectangle
(0.8 x 0.6 at z = 0.47, area 0.48).  The rectangle is named ``PortSheet`` and
the volume ``AirBox`` two ways each -- ``setEntityName`` and a named physical
group -- then the model is written to STEP, the session finalized, and the
file re-imported in a fresh session.  Every dim-2 / dim-3 entity name, bbox
and ``getMass`` is printed before and after ``occ.removeAllDuplicates()``.

``write_named_box_port_step(path)`` is the writer the step-1d test imports.
"""

import os
import sys
import tempfile

import gmsh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe_geo37_interior_facet import RECT, RECT_AREA  # noqa: E402

PORT_NAME = "PortSheet"
VOL_NAME = "AirBox"


def write_named_box_port_step(path):
    """Write the named box + rectangle STEP to ``path`` (own gmsh session).

    Must be called with no gmsh session open; finalizes its own session.
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("geo34_1d_writer")
        occ = gmsh.model.occ
        box = occ.addBox(0, 0, 0, 1, 1, 1)
        rect = occ.addRectangle(*RECT)
        occ.fragment([(3, box)], [(2, rect)])
        occ.synchronize()
        vols = gmsh.model.getEntities(3)
        if len(vols) != 1:
            raise RuntimeError(f"expected 1 volume, got {vols}")
        rect_tags = []
        for _, s in gmsh.model.getEntities(2):
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, s)
            if abs(z1 - z0) < 1e-6 and abs(0.5 * (z0 + z1) - RECT[2]) < 1e-6:
                rect_tags.append(s)
        if len(rect_tags) != 1:
            raise RuntimeError(f"expected 1 interior face, got {rect_tags}")
        s = rect_tags[0]
        if abs(occ.getMass(2, s) - RECT_AREA) > 1e-12:
            raise RuntimeError("interior face area != 0.48")
        v = vols[0][1]
        gmsh.model.setEntityName(2, s, PORT_NAME)
        gmsh.model.addPhysicalGroup(2, [s], 201)
        gmsh.model.setPhysicalName(2, 201, PORT_NAME)
        gmsh.model.setEntityName(3, v, VOL_NAME)
        gmsh.model.addPhysicalGroup(3, [v], 1)
        gmsh.model.setPhysicalName(3, 1, VOL_NAME)
        gmsh.write(path)
    finally:
        gmsh.finalize()


def _dump(label):
    print(f"--- {label} ---", flush=True)
    for dim in (3, 2):
        ents = gmsh.model.getEntities(dim)
        print(f"dim {dim}: {len(ents)} entities", flush=True)
        for _, t in ents:
            name = gmsh.model.getEntityName(dim, t)
            bb = gmsh.model.getBoundingBox(dim, t)
            mass = gmsh.model.occ.getMass(dim, t)
            print(
                f"  ({dim},{t}) name={name!r} mass={mass!r} "
                f"bbox={tuple(round(x, 9) for x in bb)}",
                flush=True,
            )
    for dim in (2, 3):
        pgs = gmsh.model.getPhysicalGroups(dim)
        print(
            f"physical groups dim {dim}: "
            f"{[(d, t, gmsh.model.getPhysicalName(d, t)) for d, t in pgs]}",
            flush=True,
        )
    vols = gmsh.model.getEntities(3)
    for _, v in vols:
        bnd = gmsh.model.getBoundary([(3, v)], oriented=False, combined=False)
        print(f"  volume {v} boundary faces: {sorted(t for _, t in bnd)}", flush=True)
        try:
            _, emb = None, gmsh.model.mesh.getEmbedded(3, v)
            print(f"  volume {v} embedded: {emb}", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"  getEmbedded failed: {exc}", flush=True)
    for _, s in gmsh.model.getEntities(2):
        if abs(gmsh.model.occ.getMass(2, s) - RECT_AREA) < 1e-12:
            up, _ = gmsh.model.getAdjacencies(2, s)
            print(f"  0.48 face ({2},{s}) upward adjacencies (volumes): " f"{list(up)}", flush=True)


def main():
    d = tempfile.mkdtemp(prefix="geo34_1d_")
    path = os.path.join(d, "box_port.step")
    write_named_box_port_step(path)
    print(f"wrote {path} ({os.path.getsize(path)} bytes)", flush=True)
    with open(path) as fh:
        txt = fh.read()
    print(f"'PortSheet' occurrences in STEP text: {txt.count(PORT_NAME)}", flush=True)
    print(f"'AirBox' occurrences in STEP text: {txt.count(VOL_NAME)}", flush=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("geo34_1d_reader")
        gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        _dump("after importShapes")
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        _dump("after removeAllDuplicates")
    finally:
        gmsh.finalize()
    print(f"gmsh version {gmsh.__version__}", flush=True)


if __name__ == "__main__":
    main()
