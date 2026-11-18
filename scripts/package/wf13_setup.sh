#!/usr/bin/env bash
# WF-13 step 1 (PROJECT_PLAN §9 item 220, brief docs/planning/chunks/WF-13.md):
# install fem-em-solver with `pip install .` into a fresh venv OUTSIDE
# /workspace, so the closed-form gates can be re-run against the installed
# package with no PYTHONPATH and a cwd outside the repository.
#
# Runs inside the fem-em-solver container; writes only under /tmp/wf13.
#   real    (first call): fresh venv + stack .pth + install + test copies
#   complex (second call, after `source /usr/local/bin/dolfinx-complex-mode`):
#           regenerate the stack .pth only (the install is pure Python).
#
# The stack .pth lists the image interpreter's own sys.path under the current
# mode (cwd=/tmp), minus '' , anything under /workspace, and the other build's
# dolfinx entries. Plain path lines: .pth files inside those directories
# (incl. the image's editable-install hook in /dolfinx-env) are NOT processed,
# so /workspace/src never re-enters.
set -euo pipefail

MODE="${1:-}"
ROOT=/tmp/wf13
VENV="$ROOT/venv"
PKG="$ROOT/pkg"
RUN="$ROOT/run"
IMAGE_PY=python3   # the image interpreter (/dolfinx-env)

case "$MODE" in
  real) OTHER=dolfinx-complex ;;
  complex) OTHER=dolfinx-real ;;
  *) echo "usage: wf13_setup.sh <real|complex>" >&2; exit 2 ;;
esac

t0=$(date +%s.%N)
step() { # step <label> <start>
  local now; now=$(date +%s.%N)
  awk -v a="$2" -v b="$now" -v l="$1" 'BEGIN { printf "[wf13] step %-28s %7.2f s\n", l, b - a }'
}

write_pth() {
  local site pth
  site=$("$VENV/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
  pth="$site/_wf13_image_stack.pth"
  ( cd /tmp && "$IMAGE_PY" - "$OTHER" > "$pth" <<'PY'
import sys
other = sys.argv[1]
for p in sys.path:
    if p == "" or p.startswith("/workspace") or other in p:
        continue
    print(p)
PY
  )
  echo "[wf13] stack .pth ($MODE): $pth"
  sed 's/^/[wf13]   /' "$pth"
}

if [ "$MODE" = real ]; then
  s=$(date +%s.%N)
  rm -rf "$ROOT"
  mkdir -p "$PKG" "$RUN/tests/validation" "$RUN/tests/package"
  echo "[wf13] (c2) before, ordinary environment (cwd=/workspace, image PYTHONPATH):"
  ( cd /workspace && "$IMAGE_PY" -c 'import sys, fem_em_solver; print("[wf13]   sys.prefix =", sys.prefix); print("[wf13]   sys.base_prefix =", sys.base_prefix); print("[wf13]   fem_em_solver.__file__ =", fem_em_solver.__file__)' )
  step "before-reading" "$s"

  s=$(date +%s.%N)
  cp /workspace/pyproject.toml /workspace/README.md "$PKG/"
  tar -C /workspace --exclude='__pycache__' --exclude='*.egg-info' -cf - src | tar -C "$PKG" -xf -
  step "copy-source" "$s"

  s=$(date +%s.%N)
  "$IMAGE_PY" -m venv --system-site-packages "$VENV"
  step "venv" "$s"

  s=$(date +%s.%N)
  write_pth
  step "stack-pth" "$s"

  s=$(date +%s.%N)
  if ! (cd /tmp && env -u PYTHONPATH "$VENV/bin/python" -c 'import setuptools; print("[wf13] setuptools", setuptools.__version__, "at", setuptools.__file__)'); then
    echo "[wf13] NEGATIVE RESULT 1: setuptools absent in the venv's view -- stop (BLOCKED)" >&2
    exit 3
  fi
  step "setuptools-probe" "$s"

  s=$(date +%s.%N)
  ( cd "$PKG" && env -u PYTHONPATH "$VENV/bin/python" -m pip install --no-index --no-deps --no-build-isolation . )
  step "pip-install" "$s"

  s=$(date +%s.%N)
  echo "[wf13] pip check (printed, not asserted):"
  ( cd /tmp && env -u PYTHONPATH "$VENV/bin/python" -m pip check ) || echo "[wf13] pip check exit status $? (ignored)"
  step "pip-check" "$s"

  s=$(date +%s.%N)
  for rel in tests/complex_mode.py \
             tests/validation/test_lossy_plane_wave.py \
             tests/validation/test_helmholtz_magnitude.py \
             tests/package/test_wf13_installed_location.py; do
    cp "/workspace/$rel" "$RUN/$rel"
    sha256sum "/workspace/$rel" "$RUN/$rel" | sed 's/^/[wf13] sha256 /'
  done
  step "copy-tests" "$s"
else
  [ -x "$VENV/bin/python" ] || { echo "[wf13] no venv at $VENV; run 'real' first" >&2; exit 2; }
  s=$(date +%s.%N)
  write_pth
  step "stack-pth" "$s"
fi

step "total-$MODE" "$t0"
