"""`WF-17` step 34 (§9 item 887, §10 U136): the sphere transient page's section
"The same two runs from YAML files" as a test.

Brief: docs/planning/chunks/WF-17.md "§9 item 887".  The solve node runs in the
**complex** build at ``-n 2``, only with ``FEM_EM_WF17_S34=1``; the other nodes
build no mesh.

Gates (bands and records imported through item 873's module, none set here):
(X)   exactly one ``sphere-transient-yaml`` block, ``run_study_file(`` once, no
      ``Study(``, no ``probe``; ``FILES`` equals item 873's two names and both
      exist under ``tests/package/wf1/``; no ```` ```yaml ```` fence on the page;
(Ref) the ``sphere-transient-yaml-check`` file: ``load_study`` raises at
      ``thermal.duty``; ``cli.main(["check", …])`` rc 2 on every rank;
(D)   ``docrefs --docs-root docs/guide --check-links`` dead 0, ``MIN_GUIDE_PAGES``;
(R)   W1: the block exec'd from ``logs/wf17/s34/`` — A run (A_t) <= ``A_BAND``
      and ``.4e`` == ``MAT12_A_T_PRINT``; S run (S2) <= ``S2_BAND`` and ``.4e``
      == ``MAT12_S2_PRINT``; sub cells == ``MAT10_SUB_CELLS``; (E) <= ``E_BAND``
      from each ``summary["temperature"]``; once per run; spy counts 2.
Controls (solve-free, asserted): (c-page) a removed or doubled tag raises
``LookupError``; (c-F) a planted ``FILES`` without the S file fails (X)'s
file-set check.

Not claimed: degree 2, duty / perfusion / probes on the sphere, any band.
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import time

import pytest
from mpi4py import MPI

import fem_em_solver
from fem_em_solver import StudyConfigError, cli, load_study
from tests.complex_mode import complex_only
from tests.package.test_wf1_yaml_sphere_transient import (
    A_BAND,
    A_FILE,
    DUTY_FILE,
    E_BAND,
    MAT10_SUB_CELLS,
    MAT12_A_T_PRINT,
    MAT12_S2_PRINT,
    S2_BAND,
    S_FILE,
    _a_t,
    _s2,
)
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf17_sphere_transient_tutorial import MIN_GUIDE_PAGES, PAGE
from tests.package.test_wf17_thermal_box_wall_tutorial import _files_of
from tests.unit.test_doc_reference_exit_codes import run_checker

REPO = PAGE.parents[2]
PAGE_ID = "docs/guide/sphere_heating_transient.md"
WF1_DIR = A_FILE.parent
SCRATCH = REPO / "logs" / "wf17" / "s34"
S34_ENV = "FEM_EM_WF17_S34"
TAG = "sphere-transient-yaml"
CHECK_TAG = "sphere-transient-yaml-check"
WANT_FILES = [A_FILE.name, S_FILE.name]


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S34_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S34_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _text() -> str:
    return PAGE.read_text()


def _file_set_ok(files) -> bool:
    """(X)'s file-set check: item 873's two names, in order, each shipped."""
    return list(files) == WANT_FILES and all((WF1_DIR / n).is_file() for n in files)


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s34_block_tags_and_files():
    text = _text()
    block = extract_tagged_block(text, TAG)
    files = _files_of(block)
    _p0(
        f"[WF-17 s34] (X) {TAG}: run_study_file( count {block.count('run_study_file(')}; "
        f"Study( {'Study(' in block.replace('run_study_file(', '')}; probe {'probe' in block}; "
        f"FILES {files}; yaml fences {text.count('```yaml')}"
    )
    assert block.count("run_study_file(") == 1
    assert "Study(" not in block.replace("run_study_file(", "")
    assert "probe" not in block
    assert _file_set_ok(files), files
    assert "```yaml" not in text
    check = extract_tagged_block(text, CHECK_TAG, lang="bash")
    assert check.strip() == f"femem check {DUTY_FILE.name}"
    # (c-F), asserted: the S file dropped from FILES fails the file-set check
    planted = block.replace(f', "{S_FILE.name}"', "")
    assert planted != block
    planted_files = _files_of(planted)
    _p0(
        f"[WF-17 s34] (c-F) planted FILES {planted_files}: file-set ok "
        f"{_file_set_ok(planted_files)}"
    )
    assert not _file_set_ok(planted_files)
    # (c-page), asserted: a removed / doubled tag raises
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {TAG}\n", "```python other\n"), TAG)
    doubled = text + "\n" + f"```python {TAG}\n" + block + "```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    _p0("[WF-17 s34] (c-page) removed and doubled tag both raise LookupError")


def test_wf17_s34_check_block_refuses():
    comm = MPI.COMM_WORLD
    try:
        load_study(DUTY_FILE)
    except StudyConfigError as exc:
        got = exc.field
    else:
        got = None
    rcs = comm.allgather(cli.main(["check", str(DUTY_FILE)]))
    _p0(f"[WF-17 s34] (Ref) {DUTY_FILE.name}: field {got!r}; rc per rank {rcs}")
    assert got == "thermal.duty"
    assert all(r == 2 for r in rcs)


def test_wf17_s34_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s34] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} rc={rc} "
        f"(MIN_GUIDE_PAGES {MIN_GUIDE_PAGES})"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------------ W1 (solve)


@complex_only
def test_wf17_s34_block_reproduces_u131(monkeypatch):
    """W1, complex build, ``-n 2``: the block's A / S runs against U131's records."""
    if not _knob():
        pytest.skip(f"{S34_ENV} unset: the sphere transient YAML block solve is opt-in")
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()

    # exec the page block from the scratch (item 838's form)
    block = extract_tagged_block(PAGE.read_bytes().decode(), TAG)
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
    t_block = time.perf_counter() - tb
    _p0(
        f"[WF-17 s34] block stdout (rank 0): {buf.getvalue()!r}; spy calls {calls['n']}; "
        f"{t_block:.1f} s"
    )
    runs = ns_["runs"]
    runA, runS = runs[A_FILE.name], runs[S_FILE.name]

    once = {k: r.study.result is r.result for k, r in runs.items()}
    n_sub = {k: int(r.thermal.mesh.topology.index_map(3).size_global) for k, r in runs.items()}
    cells = {k: r.summary["cells"] for k, r in runs.items()}
    err_A, q_cf = _a_t(runA.thermal)
    s2, ratio = _s2(runS.thermal, comm)
    e = {k: r.summary["temperature"]["max_balance_residual"] for k, r in runs.items()}
    _p0(
        f"[WF-17 s34] (Once) {once}; spy {calls['n']}; sub cells {n_sub} "
        f"(MAT10_SUB_CELLS {MAT10_SUB_CELLS}); summary cells {cells}"
    )
    _p0(
        f"[WF-17 s34] (A_t) {err_A:.10e} [{err_A:.4e}] (A_BAND {A_BAND}; record "
        f"{MAT12_A_T_PRINT}); q_cf {q_cf:.6f} W/m^3"
    )
    _p0(
        f"[WF-17 s34] (S2) ratio {ratio:.9f} |ratio - 1| {s2:.10e} [{s2:.4e}] "
        f"(S2_BAND {S2_BAND:.0e}; record {MAT12_S2_PRINT})"
    )
    _p0(f"[WF-17 s34] (E) {e} (E_BAND {E_BAND:.0e})")
    _p0(f"[WF-17 s34] elapsed {time.perf_counter() - t0:.1f} s")

    gates = {
        "spy": calls["n"] == 2 and len(runs) == 2,
        "once": all(once.values()),
        "sub_cells": all(n == MAT10_SUB_CELLS for n in n_sub.values()),
        "A_t": err_A <= A_BAND,
        "R_A_t": f"{err_A:.4e}" == MAT12_A_T_PRINT,
        "S2": s2 <= S2_BAND,
        "R_S2": f"{s2:.4e}" == MAT12_S2_PRINT,
        "E": all(v <= E_BAND for v in e.values()),
    }
    _p0(f"[WF-17 s34] gates {gates}")
    assert all(gates.values()), gates
