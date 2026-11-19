"""`GEO-34` step 1d probe (§9 item 333, chain U12a): can Gmsh's own OCC kernel
carry an interior face through a STEP file in any export form?

Measurement only (§3): asserts nothing, imported by nothing, never collected
(``python_files = test_*.py``).  One rank, plain ``python3``.

Fixture: `GEO-37` item 222's unit box and the 0.8 x 0.6 rectangle at z = 0.47
(``RECT`` / ``RECT_AREA`` from ``probe_geo37_interior_facet``).  Export forms:

- ``F-frag``  -- W0's writer ``write_named_box_port_step`` unchanged (control);
- ``F-solid`` -- the unfragmented unit box only (control);
- ``F-comp``  -- unfragmented box + free ``addRectangle(*RECT)``, one write;
- ``F-shell`` -- the rectangle alone;
- ``F-two``   -- F-solid's file and F-shell's file both imported in one reader
  session (beyond §10's two forms, review ruling Q5).

Each file is read twice, each read in its own gmsh session: r1 =
``importShapes(path)`` (the importer's read, ``step_import.py:164``), r2 =
``importShapes(path, highestDimOnly=False)``.  Counts are printed before and
after ``occ.removeAllDuplicates()``.  No ``Geometry.OCC*`` option is set.

The driver (no ``--form``) runs one child process per form and prints the
table; it exits with the largest child return code.  STEP scratch lives under
``/workspace/logs/geo34_step1d_probe/<form>/`` (gitignored).
"""

import argparse
import inspect
import os
import re
import subprocess
import sys
import time

import gmsh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe_geo34_step_names import (  # noqa: E402
    PORT_NAME,
    VOL_NAME,
    write_named_box_port_step,
)
from probe_geo37_interior_facet import RECT, RECT_AREA  # noqa: E402

FORMS = ("F-frag", "F-solid", "F-comp", "F-shell", "F-two")
SCRATCH = "/workspace/logs/geo34_step1d_probe"
OCC_OPTIONS = (
    "OCCImportLabels",
    "OCCSewFaces",
    "OCCMakeSolids",
    "OCCFixDegenerated",
    "OCCFixSmallEdges",
    "OCCFixSmallFaces",
)
RECT_BBOX = (RECT[0], RECT[1], RECT[2], RECT[0] + RECT[3], RECT[1] + RECT[4], RECT[2])


def _p(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- writers ---
def _name(dim, tag, name, pg_tag):
    gmsh.model.setEntityName(dim, tag, name)
    gmsh.model.addPhysicalGroup(dim, [tag], pg_tag)
    gmsh.model.setPhysicalName(dim, pg_tag, name)


def _write(path, with_box, with_rect):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("geo34_1d_probe_writer")
        occ = gmsh.model.occ
        box = occ.addBox(0, 0, 0, 1, 1, 1) if with_box else None
        rect = occ.addRectangle(*RECT) if with_rect else None
        occ.synchronize()
        if box is not None:
            _name(3, box, VOL_NAME, 1)
        if rect is not None:
            _name(2, rect, PORT_NAME, 201)
        gmsh.write(path)
    finally:
        gmsh.finalize()


def _write_form(form):
    d = os.path.join(SCRATCH, form)
    os.makedirs(d, exist_ok=True)
    if form == "F-frag":
        p = os.path.join(d, "box_port_frag.step")
        write_named_box_port_step(p)
        return [p]
    if form == "F-solid":
        p = os.path.join(d, "box_solid.step")
        _write(p, True, False)
        return [p]
    if form == "F-comp":
        p = os.path.join(d, "box_plus_free_rect.step")
        _write(p, True, True)
        return [p]
    if form == "F-shell":
        p = os.path.join(d, "rect_shell.step")
        _write(p, False, True)
        return [p]
    if form == "F-two":
        p1 = os.path.join(d, "box_solid.step")
        p2 = os.path.join(d, "rect_shell.step")
        _write(p1, True, False)
        _write(p2, False, True)
        return [p1, p2]
    raise ValueError(form)


# ---------------------------------------------------------------- readers ---
def _is_rect(s):
    a = gmsh.model.occ.getMass(2, s)
    bb = gmsh.model.getBoundingBox(2, s)
    return abs(a / RECT_AREA - 1.0) <= 1e-12 and all(
        abs(b - r) <= 1e-6 for b, r in zip(bb, RECT_BBOX)
    )


def _census(label):
    vols = [t for _, t in gmsh.model.getEntities(3)]
    faces = [t for _, t in gmsh.model.getEntities(2)]
    _p(f"--- {label} ---")
    _p(f"dim 3: {len(vols)} entities")
    for v in vols:
        _p(
            f"  (3,{v}) name={gmsh.model.getEntityName(3, v)!r} "
            f"mass={gmsh.model.occ.getMass(3, v)!r} "
            f"bbox={tuple(round(x, 9) for x in gmsh.model.getBoundingBox(3, v))}"
        )
    _p(f"dim 2: {len(faces)} entities")
    n_rect = n_free = n_emb = 0
    for s in faces:
        up, _ = gmsh.model.getAdjacencies(2, s)
        up = [int(u) for u in up]
        rect = _is_rect(s)
        _p(
            f"  (2,{s}) name={gmsh.model.getEntityName(2, s)!r} "
            f"mass={gmsh.model.occ.getMass(2, s)!r} "
            f"bbox={tuple(round(x, 9) for x in gmsh.model.getBoundingBox(2, s))} "
            f"up={up} rect={rect}"
        )
        if not rect:
            continue
        n_rect += 1
        if not up:
            n_free += 1
        elif len(vols) == 1 and up == vols:
            n_emb += 1
    for v in vols:
        try:
            emb = gmsh.model.mesh.getEmbedded(3, v)
            _p(f"  volume {v} getEmbedded: {[tuple(e) for e in emb]}")
        except Exception as exc:  # noqa: BLE001
            _p(f"  volume {v} getEmbedded failed: {exc}")
    for dim in (2, 3):
        pgs = gmsh.model.getPhysicalGroups(dim)
        _p(
            f"  physical groups dim {dim}: "
            f"{[(d, t, gmsh.model.getPhysicalName(d, t)) for d, t in pgs]}"
        )
    _p(
        f"counts {label}: dim3={len(vols)} dim2={len(faces)} "
        f"rect048={n_rect} free={n_free} embedded={n_emb}"
    )
    return len(vols), len(faces), n_rect, n_free, n_emb


def _read(form, paths, read):
    kwargs = {} if read == "r1" else {"highestDimOnly": False}
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(f"geo34_1d_probe_reader_{read}")
        if read == "r1":
            _p("Geometry.OCC* options (read, never set):")
            for o in OCC_OPTIONS:
                try:
                    _p(f"  Geometry.{o} = {gmsh.option.getNumber('Geometry.' + o)!r}")
                except Exception as exc:  # noqa: BLE001
                    _p(f"  Geometry.{o} unreadable: {exc}")
        for p in paths:
            out = gmsh.model.occ.importShapes(p, **kwargs)
            _p(f"importShapes({os.path.basename(p)!r}, {kwargs}) -> {out}")
        gmsh.model.occ.synchronize()
        rows = []
        for stage in ("pre", "post"):
            if stage == "post":
                gmsh.model.occ.removeAllDuplicates()
                gmsh.model.occ.synchronize()
            c = _census(f"{form} {read} {stage}-heal")
            rows.append((stage, c))
        return rows
    finally:
        gmsh.finalize()


def child(form):
    _p(f"=== form {form} ===")
    _p(f"gmsh version {gmsh.__version__}")
    _p(f"importShapes signature: {inspect.signature(gmsh.model.occ.importShapes)}")
    paths = _write_form(form)
    names = []
    for p in paths:
        with open(p) as fh:
            txt = fh.read()
        n_ent = len(re.findall(r"^#\d+\s*=", txt, flags=re.M))
        cp, ca = txt.count(PORT_NAME), txt.count(VOL_NAME)
        names.append((cp, ca))
        _p(
            f"file {p}: {os.path.getsize(p)} bytes, {n_ent} STEP entities, "
            f"'{PORT_NAME}' x{cp}, '{VOL_NAME}' x{ca}"
        )
    port_occ = sum(c[0] for c in names)
    vol_occ = sum(c[1] for c in names)
    for read in ("r1", "r2"):
        for stage, (d3, d2, nr, nf, ne) in _read(form, paths, read):
            _p(f"ROW|{form}|{read}|{stage}|OK|{d3}|{d2}|{nr}|{nf}|{ne}|" f"{port_occ}/{vol_occ}")


# ----------------------------------------------------------------- driver ---
def driver():
    rows, worst = [], 0
    for form in FORMS:
        t0 = time.perf_counter()
        r = subprocess.run(
            [sys.executable, os.path.abspath(__file__), "--form", form],
            capture_output=True,
            text=True,
        )
        dt = time.perf_counter() - t0
        sys.stdout.write(r.stdout)
        sys.stdout.flush()
        if r.stderr.strip():
            _p(f"[{form} stderr]\n{r.stderr.rstrip()}")
        _p(f"[{form}] child rc={r.returncode} elapsed {dt:.2f} s")
        worst = max(worst, abs(r.returncode))
        got = [ln[4:] for ln in r.stdout.splitlines() if ln.startswith("ROW|")]
        if r.returncode != 0 or len(got) != 4:
            tail = " / ".join(r.stderr.strip().splitlines()[-3:])
            rows.append(f"{form}|-|-|FAIL|-|-|-|-|-|-|{tail}")
        rows.extend(got)
    _p("")
    _p(
        "TABLE (form | read | stage | OK/FAIL | dim-3 | dim-2 | 0.48-faces | "
        "free | embedded | PortSheet/AirBox occurrences)"
    )
    _p("| form | read | stage | OK/FAIL | dim-3 | dim-2 | 0.48 | free | " "embedded | names |")
    _p("|---|---|---|---|---|---|---|---|---|---|")
    for row in rows:
        _p("| " + " | ".join(row.split("|")) + " |")
    _p(f"driver exit code {worst}")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--form", choices=FORMS)
    a = ap.parse_args()
    if a.form:
        child(a.form)
        return 0
    return driver()


if __name__ == "__main__":
    sys.exit(main())
