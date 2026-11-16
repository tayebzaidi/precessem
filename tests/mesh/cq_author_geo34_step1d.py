"""`GEO-34` step 1d / U12' (§9 item 985) — author the CadQuery STEP fixture.

Not collected (``python_files = test_*.py``).  Runs under ``cq-python``, which
clears ``PYTHONPATH`` (``docker/Dockerfile:41``), so the fixture constants of
``tests/mesh/probe_geo37_interior_facet.py:30-32`` are restated as literals:
RECT = (0.1, 0.2, 0.47, 0.8, 0.6), area 0.48, in a unit box.

The box is product ``AirBox``; the sheet is an explicit 4-corner polygon face
(never ``makePlane``, whose orientation is ambiguous), product ``PortSheet``.
Usage: ``cq-python tests/mesh/cq_author_geo34_step1d.py <out.step>``.
"""

import hashlib
import os
import sys

import cadquery as cq

X, Y0, Z, W, HH = 0.1, 0.2, 0.47, 0.8, 0.6
AREA = 0.48


def main(out):
    print("cadquery", cq.__version__, flush=True)
    box = cq.Workplane("XY").box(1, 1, 1, centered=False)
    corners = [
        cq.Vector(X, Y0, Z),
        cq.Vector(X + W, Y0, Z),
        cq.Vector(X + W, Y0 + HH, Z),
        cq.Vector(X, Y0 + HH, Z),
    ]
    face = cq.Face.makeFromWires(cq.Wire.makePolygon(corners, close=True))
    area = face.Area()
    bb = face.BoundingBox()
    print(f"[author] PortSheet area={area!r} bbox=({bb.xmin!r}, {bb.ymin!r}, "
          f"{bb.zmin!r}, {bb.xmax!r}, {bb.ymax!r}, {bb.zmax!r})", flush=True)
    assert abs(area - AREA) <= 1e-12, area
    # OCC pads bounding boxes by the shape tolerance (~1e-7): check at 1e-6.
    for got, want in ((bb.xmin, X), (bb.ymin, Y0), (bb.xmax, X + W),
                      (bb.ymax, Y0 + HH), (bb.zmin, Z), (bb.zmax, Z)):
        assert abs(got - want) <= 1e-6, (got, want)
    vol = box.val().Volume()
    print(f"[author] AirBox volume={vol!r}", flush=True)
    assert abs(vol - 1.0) <= 1e-12, vol

    asm = cq.Assembly(name="geo34_step1d")
    asm.add(box, name="AirBox")
    asm.add(cq.Workplane().add(face), name="PortSheet")
    if hasattr(asm, "export"):
        asm.export(out)
        route = "Assembly.export"
    else:
        asm.save(out, exportType="STEP")
        route = "Assembly.save (export absent; disclosed)"
    with open(out, "rb") as handle:
        data = handle.read()
    print(f"[author] route={route}", flush=True)
    print(f"[author] wrote {out} bytes={len(data)} "
          f"sha256={hashlib.sha256(data).hexdigest()}", flush=True)
    print(f"[author] size_on_disk={os.path.getsize(out)}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
