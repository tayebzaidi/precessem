#!/usr/bin/env bash
# install_cadquery_image.sh — rebuild the fem-em-solver image with CadQuery and
# swap the running container onto it, without disturbing a scheduled session.
#
# CadQuery is for AUTHORING test STEP files (GEO-34 step 1c, operator
# 2026-10-05). docker/Dockerfile installs it in an isolated venv
# (/opt/cadquery, run as `cq-python`) so the numpy / vtk the DolfinX builds are
# linked against are never touched. The solver reads STEP through gmsh.
#
# Run by the operator from the repo root, in a normal terminal (needs network
# for the build and the docker socket — neither is available to a sandboxed
# Claude session):
#     scripts/maintenance/install_cadquery_image.sh
#
# Steps:
#   1. record the solver environment in the OLD container (numpy, dolfinx, gmsh,
#      vtk, pip freeze) — the "before" side of the no-change check;
#   2. build the new image (the running container keeps the old one meanwhile,
#      so slots keep working during the build);
#   3. take the automation flock (waits for a running slot / review to finish;
#      new ones skip while it is held) and recreate ONLY `fem-em-solver`
#      (the xl / xxl / palace services are not touched);
#   4. verify in the new container:
#        (a) `cq-python` imports cadquery and writes a STEP file (a box plus an
#            embedded named rectangle sheet, in metres);
#        (b) the solver's Python reads that STEP through gmsh: 1 volume, the
#            sheet's face present, bounding box as drawn;
#        (c) the solver environment is unchanged: `pip freeze` and the
#            numpy / dolfinx / gmsh / vtk versions are identical to step 1,
#            and both DolfinX builds (real, complex) import;
#        (d) the STEP-import and interior-facet tests pass at -n 2
#            (tests/mesh/test_step_import.py, test_geo37_interior_facet.py).
#   Any failure exits non-zero and says which check; the old image is kept as
#   fem-em-solver:pre-cadquery so you can roll back with
#       docker tag fem-em-solver:pre-cadquery fem-em-solver:latest
#       docker compose -f docker/docker-compose.yml up -d --force-recreate fem-em-solver
#
# Output: logs/maintenance/cadquery_image_<UTC>.log (gitignored), plus the
# STEP sample at logs/maintenance/cq_sample_<UTC>.step. Nothing tracked is
# written, so the tree stays clean for the next slot.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
COMPOSE=(docker compose -f docker/docker-compose.yml)
SVC=fem-em-solver
TS="$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p logs/maintenance
LOG="logs/maintenance/cadquery_image_${TS}.log"
STEP_REL="logs/maintenance/cq_sample_${TS}.step"
LOCK="${FEM_EM_AUTOMATION_LOCK:-$HOME/.fem-em-automation.lock}"
exec > >(tee -a "$LOG") 2>&1

say() { printf '\n== %s  %s\n' "$(date '+%H:%M:%S')" "$*"; }
die() { printf '\nFAILED: %s\nlog: %s\n' "$*" "$LOG"; exit 1; }

ENV_PROBE='import importlib, json, subprocess, sys
out = {}
for m in ("numpy", "gmsh", "vtk", "mpi4py", "h5py"):
    try:
        mod = importlib.import_module(m); out[m] = getattr(mod, "__version__", None) or getattr(mod, "VTK_VERSION", "?")
    except Exception as e:
        out[m] = f"IMPORT FAILED: {e}"
out["pip_freeze"] = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout
print(json.dumps(out))'

say "1. solver environment in the current container (before)"
"${COMPOSE[@]}" ps --format '{{.Name}} {{.Image}} {{.Status}}'
BEFORE="$("${COMPOSE[@]}" exec -T "$SVC" python3 -c "$ENV_PROBE")" || die "cannot exec into $SVC (is it Up?)"
python3 -c "import json,sys; d=json.loads(sys.argv[1]); print({k:v for k,v in d.items() if k!='pip_freeze'}); print(len(d['pip_freeze'].splitlines()),'pip packages')" "$BEFORE"

say "2. build (old container keeps running)"
if docker image inspect fem-em-solver:latest >/dev/null 2>&1; then
  docker tag fem-em-solver:latest fem-em-solver:pre-cadquery
  echo "old image tagged fem-em-solver:pre-cadquery (rollback target)"
fi
"${COMPOSE[@]}" build "$SVC" || die "image build"

say "3. waiting for the automation lock ($LOCK), then recreating $SVC only"
exec 9>"$LOCK"
flock 9
echo "lock held — no scheduled session can start until this script exits"
"${COMPOSE[@]}" up -d --force-recreate --no-deps "$SVC" || die "recreate"
"${COMPOSE[@]}" ps --format '{{.Name}} {{.Image}} {{.Status}}'

say "4a. cq-python writes a STEP file"
"${COMPOSE[@]}" exec -T "$SVC" cq-python - "$STEP_REL" <<'PY' || die "4a cadquery STEP export"
import sys
import cadquery as cq
out = "/workspace/" + sys.argv[1]
print("cadquery", cq.__version__)
box = cq.Workplane("XY").box(0.10, 0.08, 0.06)                     # metres, centred at 0
sheet =cq.Face.makePlane(0.007, 0.008, basePnt=(0, 0, 0), dir=(1, 0, 0))   # port-sized rectangle, normal +x
asm = cq.Assembly(name="cq_sample")
asm.add(box, name="AirBox")
asm.add(cq.Workplane().add(sheet), name="PortSheet")
asm.save(out, exportType="STEP")
print("wrote", out)
PY

say "4b. the solver side reads it through gmsh"
"${COMPOSE[@]}" exec -T "$SVC" python3 - "$STEP_REL" <<'PY' || die "4b gmsh read of the CadQuery STEP"
import sys
import gmsh
path = "/workspace/" + sys.argv[1]
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)  # gmsh 4.15 has no -noenv
# highestDimOnly=False: the default drops the free 2D sheet beside the volume
gmsh.model.occ.importShapes(path, highestDimOnly=False); gmsh.model.occ.synchronize()
vols = gmsh.model.getEntities(3); faces = gmsh.model.getEntities(2)
bb = gmsh.model.getBoundingBox(-1, -1)
print("volumes", len(vols), "faces", len(faces), "bbox", [round(v, 6) for v in bb])
names = [gmsh.model.getEntityName(d, t) for d, t in vols + faces]
print("entity names carried by the STEP:", sorted(set(n for n in names if n)) or "(none)")
assert len(vols) == 1, f"expected 1 volume, got {len(vols)}"
assert len(faces) >= 7, f"expected the box's 6 faces plus the sheet, got {len(faces)}"
want = [-0.05, -0.04, -0.03, 0.05, 0.04, 0.03]
assert max(abs(a - b) for a, b in zip(bb, want)) < 1e-6, f"bbox {bb} != {want} (units?)"
gmsh.finalize()
print("gmsh read OK")
PY

say "4c. solver environment unchanged; both DolfinX builds import"
AFTER="$("${COMPOSE[@]}" exec -T "$SVC" python3 -c "$ENV_PROBE")" || die "4c env probe"
python3 - "$BEFORE" "$AFTER" <<'PY' || die "4c the solver environment changed (see diff above)"
import difflib, json, sys
b, a = json.loads(sys.argv[1]), json.loads(sys.argv[2])
bad = [k for k in b if k != "pip_freeze" and b[k] != a.get(k)]
for k in bad:
    print(f"  {k}: {b[k]} -> {a.get(k)}")
diff = list(difflib.unified_diff(b["pip_freeze"].splitlines(), a["pip_freeze"].splitlines(), lineterm="", n=0))
print("\n".join(diff) if diff else "pip freeze identical")
sys.exit(1 if bad or diff else 0)
PY
"${COMPOSE[@]}" exec -T "$SVC" bash -lc 'source /usr/local/bin/dolfinx-real-mode && python3 -c "import dolfinx; print(\"real\", dolfinx.__version__)"' || die "4c real build"
"${COMPOSE[@]}" exec -T "$SVC" bash -lc 'source /usr/local/bin/dolfinx-complex-mode && python3 -c "import dolfinx, numpy; from petsc4py import PETSc; assert numpy.dtype(PETSc.ScalarType).kind == \"c\"; print(\"complex\", dolfinx.__version__)"' || die "4c complex build"

say "4d. STEP-import and interior-facet tests at -n 2"
"${COMPOSE[@]}" exec -T "$SVC" bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/mesh/test_step_import.py tests/mesh/test_geo37_interior_facet.py -q --tb=short' || die "4d tests"

say "DONE — CadQuery is in the image (cq-python), solver environment unchanged"
echo "log: $LOG"
echo "sample STEP: $STEP_REL"
echo "rollback image kept as fem-em-solver:pre-cadquery (docker image rm it once you are happy)"
