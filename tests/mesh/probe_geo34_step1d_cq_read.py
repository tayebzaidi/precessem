"""`GEO-34` step 1d / U12' (§9 item 985) W0b — read the CadQuery STEP in gmsh.

Not collected.  One rank, one gmsh session, two reads (``gmsh.clear()``
between): r1 ``importShapes`` with defaults (``highestDimOnly=True``), r2 with
``highestDimOnly=False``.  For each it prints per-entity names and dims before
and after ``removeAllDuplicates``; for r2 it then fragments every volume
against the faces named ``PortSheet`` and prints the ``fragment`` map with
each child's area, bbox (printed only, never used to select) and upward
adjacency, and whether the child is in the combined wall set.
"""

import hashlib
import sys

import gmsh


def _names(tag):
    out = []
    for dim in range(4):
        for _, t in gmsh.model.getEntities(dim):
            name = gmsh.model.getEntityName(dim, t)
            out.append((dim, t, name))
            print(f"[W0b] {tag} dim={dim} tag={t} name={name!r}", flush=True)
    return out


def _named_faces(entries):
    return [
        t
        for d, t, n in entries
        if d == 2 and n and (n == "PortSheet" or n.split("/")[-1] == "PortSheet")
    ]


def read(path, hdo):
    label = f"r{'1' if hdo else '2'}(highestDimOnly={hdo})"
    gmsh.clear()
    gmsh.model.add(label)
    gmsh.model.occ.importShapes(path, highestDimOnly=hdo)
    gmsh.model.occ.synchronize()
    pre = _names(f"{label} pre-heal")
    print(f"[W0b] {label} pre-heal named PortSheet faces: {_named_faces(pre)}", flush=True)
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    post = _names(f"{label} post-heal")
    named = _named_faces(post)
    print(f"[W0b] {label} post-heal named PortSheet faces: {named}", flush=True)
    print(
        f"[W0b] {label} counts post-heal: vols={len(gmsh.model.getEntities(3))} "
        f"faces={len(gmsh.model.getEntities(2))}",
        flush=True,
    )
    if not hdo and named:
        masses = {s: gmsh.model.occ.getMass(2, s) for s in named}
        print(f"[W0b] {label} pre-fragment masses {masses!r}", flush=True)
        vols = gmsh.model.getEntities(3)
        out, out_map = gmsh.model.occ.fragment(vols, [(2, s) for s in named])
        gmsh.model.occ.synchronize()
        print(f"[W0b] {label} fragment out={out}", flush=True)
        objs = list(vols) + [(2, s) for s in named]
        for src, children in zip(objs, out_map):
            print(f"[W0b] {label} out_map {src} -> {children}", flush=True)
        walls = {
            abs(s)
            for d, s in gmsh.model.getBoundary(
                gmsh.model.getEntities(3), combined=True, oriented=False, recursive=False
            )
            if d == 2
        }
        print(f"[W0b] {label} wall set {sorted(walls)}", flush=True)
        for (_, s), children in zip([(2, s) for s in named], out_map[len(vols) :]):
            for d, c in children:
                up, _down = gmsh.model.getAdjacencies(d, c)
                print(
                    f"[W0b] {label} child of {s}: ({d},{c}) name="
                    f"{gmsh.model.getEntityName(d, c)!r} area={gmsh.model.occ.getMass(d, c)!r} "
                    f"|area-premass|={abs(gmsh.model.occ.getMass(d, c) - masses[s]):.3e} "
                    f"bbox={gmsh.model.getBoundingBox(d, c)} upward={list(up)} "
                    f"in_wall_set={c in walls}",
                    flush=True,
                )
        for _, v in gmsh.model.getEntities(3):
            print(
                f"[W0b] {label} getEmbedded(3,{v})={gmsh.model.mesh.getEmbedded(3, v)} "
                f"boundary(3,{v})={gmsh.model.getBoundary([(3, v)], oriented=False)}",
                flush=True,
            )
        # Mesh-level embedding: a conformal mesh shares the sheet's nodes with
        # the volume's tetrahedra.  Printed, not used to select.
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.08)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.08)
        gmsh.model.mesh.generate(3)
        for s in named:
            _, etags, enodes = gmsh.model.mesh.getElements(2, s)
            tri_nodes = {int(n) for n in enodes[0]} if enodes else set()
            _, _, tet_nodes = gmsh.model.mesh.getElements(3, 1)
            tet_set = {int(n) for n in tet_nodes[0]} if tet_nodes else set()
            print(
                f"[W0b] {label} mesh: sheet {s} n_tri={sum(len(t) for t in etags)} "
                f"sheet nodes={len(tri_nodes)} of which in tets={len(tri_nodes & tet_set)}",
                flush=True,
            )
        _names(f"{label} post-fragment")


def main(path):
    with open(path, "rb") as handle:
        print(f"[W0b] {path} sha256={hashlib.sha256(handle.read()).hexdigest()}", flush=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    try:
        read(path, True)
        read(path, False)
        print(f"[W0b] model bbox={gmsh.model.getBoundingBox(-1, -1)}", flush=True)
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv[1])
