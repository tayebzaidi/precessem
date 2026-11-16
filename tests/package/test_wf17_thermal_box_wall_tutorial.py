"""`WF-17` step 31 (§9 item 838, §10 U126): the thermal-box page's section
"Start from a gradient, inside an acrylic wall" as a test.

Brief: docs/planning/chunks/WF-17.md "§9 item 838" (and the review ruling
2026-11-05 03:00 above item 827).  The solve node runs in the **real** build
at ``-n 2``, only with ``FEM_EM_WF17_S31=1``; the other nodes build no mesh.

Anchors (bands imported, unmoved): the ``thermal-box-wall-gradient`` block,
exec'd from ``logs/wf17/s31/``, gives θ(P5) on the R1 file equal to P63's
direct R1 (computed in the same window from item 830's module's helpers,
imported by name) at ``IDENTITY_BAND``; (E) ≤ ``E_BAND`` on both runs;
the ``thermal-box-wall-check`` block's file returns rc 2 and names
``geometry.wall.thickness_m``; (D) 0 dead docrefs, ``MIN_GUIDE_PAGES`` holds.
Controls (asserted, rule (e)): (Drift) the R0 file's θ(P5) misses R1 by P63
W1b's printed drift, +1.828983e-01 K, within 1e-7
(``20261104T060527Z_MAT-21-P63-W1b.log:76``); (c-page) a removed or doubled tag
raises ``LookupError``.

Not claimed: anything beyond the page section (the six-wall box, a wall under
duty, a cited PMMA ρc).
"""

from __future__ import annotations

import ast
import contextlib
import io
import os
import shutil
import time
from pathlib import Path

import pytest
from mpi4py import MPI

import fem_em_solver
from fem_em_solver import cli, load_study
from fem_em_solver.thermal import solve_transient_conduction

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf17_first_study_tutorial import LINK_PASS_RE, _checker_stdout, extract_tagged_block
from tests.package.test_wf1_yaml_thermal_box_wall import (
    CELLS,
    D_WALL,
    DRIFT_BOUND,
    DT,
    G_GRAD,
    H_GEL,
    I_P5,
    IDENTITY_BAND,
    MALFORMED,
    N_MAIN,
    Q,
    R0_FILE,
    R1_FILE,
    READ_PTS,
    ROBIN,
    RUNG,
    T_INF,
    W1B_DRIFT_P5,
    A_BOX,
    B_BOX,
    C_BOX,
    E_BAND,
    _field,
    _materials,
    _read,
    _rho_c_field,
    _robin_facets,
    _walled_box,
)
from tests.unit.test_doc_reference_exit_codes import run_checker

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "thermal_box_from_yaml.md"
PAGE_ID = "docs/guide/thermal_box_from_yaml.md"
WF1_DIR = REPO / "tests" / "package" / "wf1"
SCRATCH = REPO / "logs" / "wf17" / "s31"
S31_ENV = "FEM_EM_WF17_S31"
TAG = "thermal-box-wall-gradient"
CHECK_TAG = "thermal-box-wall-check"
WANT_FILES = [R1_FILE.name, R0_FILE.name]
#: the latest page census floor (`test_wf17_estimate_human_tutorial.py:75`)
MIN_GUIDE_PAGES = 68


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S31_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S31_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _text() -> str:
    return PAGE.read_text()


def _files_of(block: str) -> list:
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "FILES" for t in node.targets
        ):
            return list(ast.literal_eval(node.value))
    raise LookupError("no FILES assignment in the block")


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s31_block_tags_and_files():
    text = _text()
    block = extract_tagged_block(text, TAG)
    files = _files_of(block)
    _p0(f"[WF-17 s31] (X) {TAG}: run_study_file( count {block.count('run_study_file(')}; "
        f"FILES {files}; yaml fences {text.count('```yaml')}")
    assert block.count("run_study_file(") == 1 and "Study(" not in block
    assert files == WANT_FILES
    for name in files:
        assert (WF1_DIR / name).is_file(), name
    assert "```yaml" not in text
    check = extract_tagged_block(text, CHECK_TAG, lang="bash")
    assert check.strip() == f"femem check {MALFORMED.name}"
    # (c-page), asserted: a removed / doubled tag raises
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {TAG}", "```python other"), TAG)
    doubled = text + "\n" + f"```python {TAG}\n" + block + "```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    _p0("[WF-17 s31] (c-page) removed and doubled tag both raise LookupError")


def test_wf17_s31_check_block_refuses():
    comm = MPI.COMM_WORLD
    got = _field(lambda: load_study(MALFORMED))
    rc = cli.main(["check", str(MALFORMED)])
    rcs = comm.allgather(rc)
    _p0(f"[WF-17 s31] (Ref) {MALFORMED.name}: field {got!r}; rc per rank {rcs}")
    assert got == "geometry.wall.thickness_m"
    assert all(r == 2 for r in rcs)


def test_wf17_s31_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s31] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} rc={rc} "
          f"(MIN_GUIDE_PAGES {MIN_GUIDE_PAGES})")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------------ W1 (solve)


def test_wf17_s31_block_reproduces_p63_r1(monkeypatch):
    """W1, real build, ``-n 2``: the block's R1 run against P63's direct R1."""
    if not _knob():
        pytest.skip(f"{S31_ENV} unset: the walled-box block solve is environment-selected")
    comm = MPI.COMM_WORLD
    assert not IS_COMPLEX, "the walled-box records are real-build"
    assert comm.size == 2, comm.size
    t0 = time.perf_counter()

    # P63's direct R1, the same lines as item 830's W1 node
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, _inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731
    direct = solve_transient_conduction(
        msh, ft, k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN,
        theta=0.5, startup_half_steps=2, dt=DT, rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad)
    th_d, ok_d = _read(direct.temperature, READ_PTS)
    cells_d = int(msh.topology.index_map(3).size_global)
    _p0(f"[WF-17 s31] direct R1: cells {cells_d} N_w {nw}; theta(P5) {th_d[I_P5]:+.12e} K; "
        f"(E) {direct.max_balance_residual:.3e}; {time.perf_counter() - t0:.1f} s")

    # exec the page block from the scratch
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), TAG)
    files = _files_of(block)
    if comm.rank == 0:
        if SCRATCH.exists():
            shutil.rmtree(SCRATCH)
        SCRATCH.mkdir(parents=True)
        for name in files:
            (SCRATCH / name).write_bytes((WF1_DIR / name).read_bytes())
    comm.Barrier()
    monkeypatch.chdir(SCRATCH)
    calls = {"n": 0}
    original = fem_em_solver.run_study_file

    def spy(path, *args, **kwargs):
        calls["n"] += 1
        return original(path, *args, **kwargs)

    monkeypatch.setattr(fem_em_solver, "run_study_file", spy)
    ns_ = {"__name__": "__main__"}
    buf = io.StringIO()
    tb = time.perf_counter()
    with contextlib.redirect_stdout(buf):
        exec(compile(block, f"{PAGE_ID}#{TAG}", "exec"), ns_)
    _p0(f"[WF-17 s31] block stdout (rank 0): {buf.getvalue()!r}; spy calls {calls['n']}; "
        f"{time.perf_counter() - tb:.1f} s")
    runs = ns_["runs"]
    r1, r0 = runs[R1_FILE.name], runs[R0_FILE.name]
    th1, ok1 = _read(r1.thermal.temperature, READ_PTS)
    th0, ok0 = _read(r0.thermal.temperature, READ_PTS)
    ref = th_d[I_P5]
    rel = abs(th1[I_P5] - ref) / abs(ref)
    drift = ref - th0[I_P5]
    e = {k: v.summary["temperature"]["max_balance_residual"] for k, v in runs.items()}
    _p0(f"[WF-17 s31] (I) block R1 theta(P5) {th1[I_P5]:+.12e} vs direct rel {rel:.3e} "
        f"(IDENTITY_BAND {IDENTITY_BAND:.0e}); (E) {e} (E_BAND {E_BAND:.0e}); (Drift) "
        f"R1 - R0 at P5 {drift:+.9e} K vs W1b {W1B_DRIFT_P5:+.6e}, |delta| "
        f"{abs(drift - W1B_DRIFT_P5):.3e} (bound {DRIFT_BOUND:g})")
    gates = {
        "valid_points": all(ok_d) and all(ok1) and all(ok0),
        "runs": len(runs) == 2 and calls["n"] == 2,
        "Cells": cells_d == CELLS and all(v.summary["cells"] == CELLS for v in runs.values()),
        "I": rel <= IDENTITY_BAND,
        "E_direct": float(direct.max_balance_residual) <= E_BAND,
        "E_runs": all(v <= E_BAND for v in e.values()),
        "Drift": abs(drift - W1B_DRIFT_P5) <= DRIFT_BOUND,
    }
    _p0(f"[WF-17 s31] gates {gates}; elapsed {time.perf_counter() - t0:.1f} s")
    assert all(gates.values()), gates
