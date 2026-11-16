"""`WF-17` step 11: the duty-cycle / container-wall tutorial page's two blocks
reproduce the `WF-15` step 3 transient records through the Study API.

Brief: docs/planning/chunks/WF-17.md, "§9 item 605".  Real build, ``-n 2``.
The duty block is `MAT-17` case A (``DutyCycle``), the wall block `MAT-18`
step 2(b) (T1) (``ContainerWall``).  The solve nodes run only with
``FEM_EM_WF17_S11=duty`` / ``=wall`` (and then assert the real build and
``-n 2``, never skip); the other nodes build no mesh.

Anchors: (X) one block per tag, one ``Study(`` call, one ``.temperature(``
call, no solve call (by AST name), the ``duty`` / ``wall`` keyword; (In) the
block's inputs equal the gate constants exactly, read back under a patched
``Study.temperature`` (no solve); (UT) the block's mesh equals the gate's
bitwise; (G) interior rel L2 <= ``INTERIOR_BAND`` and balance <=
``BALANCE_BAND``; (R) the records at ``.4e``; (D) 0 dead links.  Controls
(rule (e)): (c1) duty dropped and (c2) wall ρc dropped, asserted to miss the
band (backed by `WF-15` step 3's same-fixture reads); (c-page) by
construction.  Every record and band is imported by name.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import os
import re
import resource
import shutil
import time
from pathlib import Path

import pytest
from mpi4py import MPI

from fem_em_solver import (
    ContainerWall,
    DutyCycle,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    TransientSpec,
    UniformHeatSource,
)

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf15_study_transient import (
    BALANCE_BAND,
    D_A,
    INTERIOR_BAND,
    MAT17_A_PRINT,
    MAT18_C1_UNIFORM_RHOC,
    MAT18_T1_CELLS,
    MAT18_T1_PRINT,
    P_A,
    RES,
    _duty_study,
    _exact_duty,
    _on,
    _same_mesh,
    _wall_study,
)
from tests.package.test_wf17_first_study_tutorial import _checker_stdout, extract_tagged_block
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_mat9_transient_conduction_sphere import (
    H,
    K,
    Q,
    R,
    RHO_C,
    T_INF,
    _ball,
    _interior_error,
)
from tests.validation.test_mat17_duty_cycle_fem import T_END_A, _n
from tests.validation.test_mat18_two_layer_sphere_closed_form import (
    A_,
    B_,
    K1,
    K2,
    Q1,
    RC1,
    RC2,
)
from tests.validation.test_mat18_two_layer_sphere_transient_fem import (
    CELL_GUARD,
    N_CN,
    T_STAR,
    _errors,
    _fine,
)

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "duty_cycle_and_container_wall.md"
PAGE_ID = "docs/guide/duty_cycle_and_container_wall.md"
SCRATCH = REPO / "logs" / "wf17" / "s11"
S11_ENV = "FEM_EM_WF17_S11"
MIN_GUIDE_PAGES = 60
RANKS = 2

DUTY_TAG = "duty-python"
WALL_TAG = "wall-python"
SOLVE_NAMES = {"solve", "run_study_file"}
YAMLS = ("thermal_ball_duty_a.yaml", "thermal_ball_wall.yaml")

# `WF-15` step 3's full-precision reads, printed beside ours, never asserted.
WF15_DUTY_FULL = "4.3869840018e-04"  # 20261008T213834Z_WF-15-s3-W1a-duty.log.gz:91
WF15_WALL_FULL = "8.2066325346e-04"  # 20261008T214046Z_WF-15-s3-W1b-wall.log.gz:91
C1_BACKING = "9.7087e-01"  # 20261008T213834Z_WF-15-s3-W1a-duty.log.gz:98

LINK_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _rss_gib(comm):
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(kib, op=MPI.SUM) / 1024.0**2


def _text() -> str:
    return PAGE.read_text()


def _call_name(node):
    fn = node.func
    return fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)


def _calls_named(block: str, name: str) -> list:
    return [n for n in ast.walk(ast.parse(block))
            if isinstance(n, ast.Call) and _call_name(n) == name]


def _attr_calls(block: str, name: str) -> list:
    return [n for n in ast.walk(ast.parse(block))
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
            and n.func.attr == name]


def _x_defect(block: str, tag: str):
    """None when the block has the (X) shape for ``tag``, else the reason."""
    n_study = len(_calls_named(block, "Study"))
    temps = _attr_calls(block, "temperature")
    solves = sorted({_call_name(n) for n in ast.walk(ast.parse(block))
                     if isinstance(n, ast.Call)} & SOLVE_NAMES)
    if n_study != 1:
        return f"{n_study} Study( calls; need one"
    if len(temps) != 1:
        return f"{len(temps)} .temperature( calls; need one"
    if solves:
        return f"solve-shaped calls {solves}"
    kws = {k.arg for k in temps[0].keywords}
    if "transient" not in kws:
        return f".temperature( keywords {sorted(kws)}: missing 'transient'"
    if tag == DUTY_TAG and "duty" not in kws:
        return f".temperature( keywords {sorted(kws)}: missing 'duty'"
    if tag == WALL_TAG:
        geos = _calls_named(block, "ThermalBallGeometry")
        gk = {k.arg for g in geos for k in g.keywords}
        if len(geos) != 1 or "wall" not in gk:
            return f"ThermalBallGeometry( calls {len(geos)}, keywords {sorted(gk)}: missing 'wall'"
    return None


def _scratch(comm):
    if comm.rank == 0:
        if SCRATCH.exists():
            shutil.rmtree(SCRATCH)
        SCRATCH.mkdir(parents=True)
    comm.Barrier()


def _exec_block(block: str, tag: str):
    ns = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{tag}", "exec")
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(buf):
        exec(code, ns)
    return ns, time.perf_counter() - t0, buf.getvalue()


def _knob():
    raw = os.environ.get(S11_ENV, "").strip()
    if raw not in ("", "duty", "wall"):
        raise ValueError(f"{S11_ENV}={raw!r}: only unset/empty, 'duty' or 'wall' are accepted")
    return raw


def _require_window(which: str):
    if _knob() != which:
        pytest.skip(f"{S11_ENV} != {which!r}: the tutorial {which} window is "
                    "environment-selected (real build, -n 2)")
    comm = MPI.COMM_WORLD
    assert not IS_COMPLEX, "the records are real-build digits; run in the real build"
    assert comm.size == RANKS, comm.size
    return comm


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s11_blocks_and_inputs(monkeypatch, tmp_path):
    comm = MPI.COMM_WORLD
    text = _text()
    blocks = {tag: extract_tagged_block(text, tag) for tag in (DUTY_TAG, WALL_TAG)}

    # (X)
    for tag, block in blocks.items():
        reason = _x_defect(block, tag)
        _p0(f"[WF-17 s11] (X) {tag}: Study( {len(_calls_named(block, 'Study'))}, "
            f".temperature( {len(_attr_calls(block, 'temperature'))}; defect {reason}")
        assert reason is None, reason

    # (In): the inputs read back under a patched Study.temperature (no solve)
    calls = []

    def _record(self, cooling, *, transient=None, duty=None, **kw):
        calls.append(dict(study=self, cooling=cooling, transient=transient, duty=duty, kw=kw))
        return None

    _scratch(comm)
    monkeypatch.chdir(SCRATCH)
    monkeypatch.setattr(Study, "temperature", _record)
    got = {}
    for tag, block in blocks.items():
        n0 = len(calls)
        ns, _, out = _exec_block(block, tag)
        assert len(calls) == n0 + 1, f"{tag}: {len(calls) - n0} temperature calls"
        c = calls[-1]
        assert ns["study"] is c["study"] and ns["result"] is None
        assert c["kw"] == {}, c["kw"]
        got[tag] = c
        _p0(f"[WF-17 s11] (In) {tag}: geometry {c['study'].geometry!r}; source "
            f"{c['study'].excitation!r}; cooling {c['cooling']!r}; transient "
            f"{c['transient']!r}; duty {c['duty']!r}; stdout {out!r}")
    monkeypatch.undo()

    d = got[DUTY_TAG]
    g, s, co, tr, du = (d["study"].geometry, d["study"].excitation, d["cooling"],
                        d["transient"], d["duty"])
    assert type(g) is ThermalBallGeometry and g.radius == R and g.resolution == RES
    assert g.wall is None
    assert type(s) is UniformHeatSource and s.q_w_per_m3 == Q
    assert type(co) is RobinCooling
    assert co.k_w_per_m_k == K and co.h_w_per_m2_k == H and co.t_inf_c == T_INF
    assert type(tr) is TransientSpec
    assert tr.dt_s == 150.0 and tr.n_steps == _n(T_END_A, 150.0) and tr.rho_c_j_per_m3_k == RHO_C
    assert type(du) is DutyCycle and du.period_s == P_A and du.duty == D_A
    # and the whole dataclasses equal `_duty_study`'s (defaults included)
    assert (g, s, co, tr) == (ThermalBallGeometry(R, RES), UniformHeatSource(Q),
                              RobinCooling(K, H, T_INF),
                              TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C))
    assert du == DutyCycle(P_A, D_A)

    w = got[WALL_TAG]
    g, s, co, tr, du = (w["study"].geometry, w["study"].excitation, w["cooling"],
                        w["transient"], w["duty"])
    assert type(g) is ThermalBallGeometry and g.radius == A_ and g.resolution == RES
    assert type(g.wall) is ContainerWall
    assert g.wall.thickness_m == B_ - A_ and g.wall.k_w_per_m_k == K2
    assert g.wall.rho_c_j_per_m3_k == RC2
    assert type(s) is UniformHeatSource and s.q_w_per_m3 == Q1
    assert type(co) is RobinCooling
    assert co.k_w_per_m_k == K1 and co.h_w_per_m2_k == H and co.t_inf_c == T_INF
    assert type(tr) is TransientSpec
    assert tr.dt_s == T_STAR / N_CN and tr.n_steps == N_CN and tr.rho_c_j_per_m3_k == RC1
    assert du is None
    assert (g, s, co, tr) == (ThermalBallGeometry(A_, RES, wall=ContainerWall(B_ - A_, K2, RC2)),
                              UniformHeatSource(Q1), RobinCooling(K1, H, T_INF),
                              TransientSpec(T_STAR / N_CN, N_CN, RC1))
    _p0("[WF-17 s11] (In) both blocks' inputs equal the gate constants exactly (no solve)")

    # (c-page), by construction: removed / doubled tags raise
    for tag in (DUTY_TAG, WALL_TAG):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```python {tag}\n", "```python other\n"), tag)
        doubled = text + "\n" + f"```python {tag}\n" + blocks[tag] + "```\n"
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag)
    _p0("[WF-17 s11] (c-page) removed / doubled tags raise LookupError for both blocks")

    # (c-page), by construction: a tmp copy with duty= removed fails (X), naming it
    altered = text.replace("    duty=DutyCycle(1200.0, 0.5),\n", "", 1)
    assert altered != text
    altered_page = tmp_path / PAGE.name
    altered_page.write_text(altered)
    reason = _x_defect(extract_tagged_block(altered_page.read_text(), DUTY_TAG), DUTY_TAG)
    _p0(f"[WF-17 s11] (c-page) duty= dropped -> {reason}")
    assert reason is not None and "'duty'" in reason

    # no reading, no artifact-shaped token on the page; the YAML links resolve
    for token in (".json", ".npz", ".xdmf", ".s4p", MAT17_A_PRINT, MAT18_T1_PRINT,
                  str(MAT18_T1_CELLS), "215 542", f"{MAT18_C1_UNIFORM_RHOC:.4e}", C1_BACKING,
                  WF15_DUTY_FULL, WF15_WALL_FULL, "e-04", "GiB"):
        assert token not in text, token
    for name in YAMLS:
        assert (REPO / "tests" / "package" / "wf1" / name).is_file(), name
        assert f"(../../tests/package/wf1/{name})" in text, name


def test_wf17_s11_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s11] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- W-duty (-n 2, real)


def test_wf17_s11_duty_block_reproduces_mat17_a(monkeypatch):
    comm = _require_window("duty")
    t00 = time.perf_counter()
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), DUTY_TAG)
    _p0(f"[WF-17 s11] page sha256 {_sha(raw)}; block {DUTY_TAG} sha256 {_sha(block.encode())}")
    msh, _, _ = _ball(RES)
    _scratch(comm)
    monkeypatch.chdir(SCRATCH)
    ns, t_block, out = _exec_block(block, DUTY_TAG)
    result = ns["result"]

    ca, cb, geo = _same_mesh(result.mesh, msh, comm)
    _p0(f"[WF-17 s11] (UT) duty: cells {ca} / {cb} (ref _ball(RES)); geometry equal {geo}; "
        f"{len(result.source_multipliers)} (half-)step multipliers; block stdout {out!r}")
    assert ca == cb, "(UT) cells"
    assert geo, "(UT) geometry"

    ex = _exact_duty(msh, T_END_A, P_A, D_A)
    e = _interior_error(msh, _on(msh, result.temperature), ex)
    bal = result.max_balance_residual
    _p0(f"[WF-17 s11] (G) duty: interior_rel_L2={e:.10e} (band {INTERIOR_BAND:.0e}; WF-15 read "
        f"{WF15_DUTY_FULL}); max_balance={bal:.3e} (band {BALANCE_BAND:.0e}); block {t_block:.1f}s")
    _p0(f"[WF-17 s11] (R) duty: {e:.4e} vs record MAT17_A_PRINT {MAT17_A_PRINT}")

    c1 = _duty_study(150.0, _n(T_END_A, 150.0), None)
    e_c1 = _interior_error(msh, _on(msh, c1.temperature), ex)
    _p0(f"[WF-17 s11] (c1) no DutyCycle: interior miss {e_c1:.4e} (backing {C1_BACKING}; asserted "
        f"> {INTERIOR_BAND:.0e}); source_multipliers {c1.source_multipliers!r}")
    _p0(f"[WF-17 s11] W-duty node {time.perf_counter() - t00:.1f}s; summed ru_maxrss "
        f"{_rss_gib(comm):.2f} GiB")

    assert bal <= BALANCE_BAND, f"(G) balance {bal:.3e}"
    assert e <= INTERIOR_BAND, f"(G) {e:.4e}"
    assert f"{e:.4e}" == MAT17_A_PRINT, f"(R) {e:.4e} != {MAT17_A_PRINT}"
    assert c1.source_multipliers == [], "(c1) multipliers"
    assert e_c1 > INTERIOR_BAND, f"(c1) did not miss: {e_c1:.4e}"


# ------------------------------------------------------------- W-wall (-n 2, real)


def test_wf17_s11_wall_block_reproduces_mat18_t1(monkeypatch):
    comm = _require_window("wall")
    t00 = time.perf_counter()
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), WALL_TAG)
    _p0(f"[WF-17 s11] page sha256 {_sha(raw)}; block {WALL_TAG} sha256 {_sha(block.encode())}")
    msh, ct, _, nglob, _, _, _, theta = _fine()
    _p0(f"[WF-17 s11] cell guard: {nglob} cells (guard {CELL_GUARD})")
    assert nglob <= CELL_GUARD, "cell guard; stop before the solve"
    _scratch(comm)
    monkeypatch.chdir(SCRATCH)
    ns, t_block, out = _exec_block(block, WALL_TAG)
    result = ns["result"]

    ca, cb, geo = _same_mesh(result.mesh, msh, comm)
    _p0(f"[WF-17 s11] (UT) wall: cells {ca} / {cb} (record MAT18_T1_CELLS {MAT18_T1_CELLS}); "
        f"geometry equal {geo}; block stdout {out!r}")
    assert ca == cb == MAT18_T1_CELLS, "(UT) cells"
    assert geo, "(UT) geometry"

    e_in = _errors(msh, ct, _on(msh, result.temperature), theta)[1]
    bal = result.max_balance_residual
    _p0(f"[WF-17 s11] (G) wall: interior_rel_L2={e_in:.10e} (band {INTERIOR_BAND:.0e}; WF-15 read "
        f"{WF15_WALL_FULL}); max_balance={bal:.3e} (band {BALANCE_BAND:.0e}); block {t_block:.1f}s")
    _p0(f"[WF-17 s11] (R) wall: {e_in:.4e} vs record MAT18_T1_PRINT {MAT18_T1_PRINT}")

    c2 = _wall_study(RC1)
    e_c2 = _errors(msh, ct, _on(msh, c2.temperature), theta)[1]
    _p0(f"[WF-17 s11] (c2) wall rho_c = RC1: interior miss {e_c2:.4e} (backing "
        f"MAT18_C1_UNIFORM_RHOC {MAT18_C1_UNIFORM_RHOC:.4e}; asserted > {INTERIOR_BAND:.0e})")
    _p0(f"[WF-17 s11] W-wall node {time.perf_counter() - t00:.1f}s; summed ru_maxrss "
        f"{_rss_gib(comm):.2f} GiB")

    assert bal <= BALANCE_BAND, f"(G) balance {bal:.3e}"
    assert e_in <= INTERIOR_BAND, f"(G) {e_in:.4e}"
    assert f"{e_in:.4e}" == MAT18_T1_PRINT, f"(R) {e_in:.4e} != {MAT18_T1_PRINT}"
    assert e_c2 > INTERIOR_BAND, f"(c2) did not miss: {e_c2:.4e}"
