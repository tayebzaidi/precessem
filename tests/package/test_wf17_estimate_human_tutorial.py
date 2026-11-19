"""`WF-17` step 20a: the human-scale estimate tutorial page reproduces the
`WF-13` F-human estimates at degree 1 (row FH1) and degree 2 (row FH2) through
its two API blocks and ``femem run --estimate``, with no solve.

Brief: docs/planning/chunks/WF-17.md, "§9 item 645" (standing clauses from
"§9 item 621").  The build node runs only with ``FEM_EM_WF17_S20A=1`` (and then
asserts the complex build and ``-n 2``, never skips).  The other nodes build no
mesh.

Anchors (records imported by name, none set here):
  (X)   one block per tag; each python block has one ``load_study(`` call, one
        ``.estimate(`` call, no solve / tune / sweep / ``run_study_file`` name;
        d1 loads the tracked file by name; d2 edits exactly ``  degree: 1`` to
        ``  degree: 2``; the CLI line's tokens; no ``--ranks`` with
        ``--estimate``.
  (In)  W0, the F-human build patched to raise: each block's Study is an
        F-human Study at degree 1 / 2, no result, and the build is reached.
  (FH1) unknowns == E_RECORD, cells == HOLE_CELL_RECORD, degree 1, "FH1" used,
        not_priced == ("thermal transient",).
  (FH2) unknowns == N_DOFS_D2_RECORD, cells == HOLE_CELL_RECORD, degree 2,
        "FH2" used, the same not_priced.
  (NoSolve) ``study.result is None``; ``Study.solve`` / ``solve_sweep`` /
        ``tune`` and ``run_study_file`` patched to raise, and the block's Study
        resolves the patched objects (c-solve).
  (CLI = API) both files: rc 0 on every rank, printed fields equal the API's,
        no output directory.
  (B)   one real F-human build per rank (memoised; the page never shows it).
  (D)   0 dead links over docs/guide.
Controls: (c-page) removed / doubled tags raise; the d1 block naming
``f_human_tuned_copper.yaml`` fails (X); the d2 block left at ``degree: 1``
fails (FH2) on degree (W0) and on degree, unknowns and row (build window).
Predicted seconds and peak memory are printed only.

`WF-17` step 26 (§9 item 718, U107), additive: the ``estimate-f-clinical-d2``
block, no mesh, no knob.  (FC2p) its printed line equals row FC2 at every digit,
built from the imported literals; (Ref) ``solve()`` and ``femem run`` refuse the
copy with both tiers named; (CLI = API); control: the block left at degree 1
misses FC2's unknowns by U2 / U1 (asserted, ``UNKNOWN_MISS``).
"""

from __future__ import annotations

import ast
import contextlib
import dataclasses
import hashlib
import io
import os
import re
import shutil
import time
from pathlib import Path

import pytest
from mpi4py import MPI

from fem_em_solver.study import FHumanBirdcageGeometry, ResourceEstimate, Study
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf13_cli_estimate import _diff, _parse, _patch_solves
from tests.package.test_wf13_estimate_f_clinical import FC, no_build  # noqa: F401
from tests.package.test_wf13_estimate_f_clinical_d2 import FC2_LITERALS, UNKNOWN_MISS
from tests.package.test_wf17_first_study_tutorial import _checker_stdout, extract_tagged_block
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_wf10_f_human_tuned_b1_sar import E_RECORD, N_DOFS_D2_RECORD
from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "estimate_human_scale.md"
PAGE_ID = "docs/guide/estimate_human_scale.md"
YAML = REPO / "tests" / "package" / "wf1" / "f_human_hole_copper_tune.yaml"
D2_NAME = "f_human_hole_copper_tune_d2.yaml"
SCRATCH = REPO / "logs" / "wf17" / "s20a"
KNOB = "FEM_EM_WF17_S20A"
MIN_GUIDE_PAGES = 68
RANKS = 2
NOT_PRICED = ("thermal transient",)

D1_TAG = "estimate-f-human-d1"
D2_TAG = "estimate-f-human-d2"
CLI_TAG = "estimate-f-human-cli"
FC_TAG = "estimate-f-clinical-d2"
FC_D2_NAME = "f_clinical_15t_d2.yaml"
SOLVE_NAMES = {"solve", "solve_sweep", "tune", "run_study_file"}
D1_LINE = '"  degree: 1"'
D2_LINE = '"  degree: 2\\n"'
D2_LEFT_AT_1 = '"  degree: 1\\n"'

LINK_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")


class _BuildReached(Exception):
    pass


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _field_names() -> list:
    return [f.name for f in dataclasses.fields(ResourceEstimate)]


def _calls_named(block: str, name: str) -> list:
    out = []
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Call):
            fn = node.func
            fname = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)
            if fname == name:
                out.append(node)
    return out


def _names_used(block: str) -> set:
    out = set()
    for n in ast.walk(ast.parse(block)):
        if isinstance(n, ast.Name):
            out.add(n.id)
        elif isinstance(n, ast.Attribute):
            out.add(n.attr)
        elif isinstance(n, (ast.Import, ast.ImportFrom)):
            out.update(a.name for a in n.names)
    return out


def _x_defect(block: str, tag: str):
    """None when the python block passes (X), else the reason."""
    n_load = len(_calls_named(block, "load_study"))
    n_est = len(_calls_named(block, "estimate"))
    bad = sorted(_names_used(block) & SOLVE_NAMES)
    if n_load != 1 or n_est != 1 or bad:
        return f"{tag}: load_study {n_load}, estimate {n_est}, solve names {bad}"
    consts = {
        n.value
        for n in ast.walk(ast.parse(block))
        if isinstance(n, ast.Constant) and isinstance(n.value, str)
    }
    if tag == D1_TAG:
        arg = _calls_named(block, "load_study")[0].args[0]
        name = ast.literal_eval(arg) if isinstance(arg, ast.Constant) else None
        if name != YAML.name:
            return f"{tag}: load_study file {name!r} != {YAML.name!r}"
    else:
        if YAML.name not in consts or D2_NAME not in consts:
            return f"{tag}: reads/writes {sorted(consts)}; need {YAML.name!r} and {D2_NAME!r}"
        if D1_LINE not in block or D2_LINE not in block:
            return f"{tag}: the edit {D1_LINE} -> {D2_LINE} is not in the block"
    return None


def _knob() -> bool:
    raw = os.environ.get(KNOB, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{KNOB}={raw!r}: only unset/empty or '1' are accepted")


def _scratch(monkeypatch):
    """Rank 0 clears the scratch and copies the tracked file in; every rank chdirs."""
    comm = MPI.COMM_WORLD
    if comm.rank == 0:
        if SCRATCH.exists():
            shutil.rmtree(SCRATCH)
        SCRATCH.mkdir(parents=True)
        data = YAML.read_bytes()
        (SCRATCH / YAML.name).write_bytes(data)
        print(
            f"[WF-17 s20a] scratch {SCRATCH.relative_to(REPO)}; copied {YAML.name} "
            f"sha256 {_sha(data)}",
            flush=True,
        )
    comm.Barrier()
    monkeypatch.chdir(SCRATCH)


def _exec(block: str, tag: str) -> dict:
    ns = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{tag}", "exec")
    try:
        with contextlib.redirect_stdout(buf):
            exec(code, ns)
    finally:
        _p0(f"[WF-17 s20a] block {tag} stdout (rank 0): {buf.getvalue()!r}")
    return ns


def _patch_all_solves(monkeypatch):
    """The gates' ``_patch_solves`` plus ``solve_sweep`` / ``tune``."""
    _patch_solves(monkeypatch)

    def _boom(*a, **k):
        raise AssertionError("solve reached")

    monkeypatch.setattr(Study, "solve_sweep", _boom)
    monkeypatch.setattr(Study, "tune", _boom)
    return Study.solve, _boom


def _resolves_patched(study, solve_obj, boom) -> bool:
    return (
        type(study).solve is solve_obj
        and type(study).solve_sweep is boom
        and type(study).tune is boom
        and isinstance(study, Study)
    )


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s20a_blocks_and_inputs(monkeypatch, tmp_path):
    """(X), (In), (c-page); the build patched to raise."""
    text = _text()
    d1 = extract_tagged_block(text, D1_TAG)
    d2 = extract_tagged_block(text, D2_TAG)
    for tag, block in ((D1_TAG, d1), (D2_TAG, d2)):
        reason = _x_defect(block, tag)
        _p0(f"[WF-17 s20a] (X) {tag}: {reason or 'ok'}")
        assert reason is None, reason
    assert YAML.is_file()
    assert sum(ln.rstrip("\n") == "  degree: 1" for ln in YAML.read_text().splitlines()) == 1
    assert f"(../../tests/package/wf1/{YAML.name})" in text

    cli = extract_tagged_block(text, CLI_TAG, lang="bash")
    lines = [ln for ln in cli.splitlines() if ln.strip()]
    _p0(f"[WF-17 s20a] (X) {CLI_TAG}: {lines!r}")
    assert len(lines) == 1
    toks = lines[0].split()
    assert toks == ["mpiexec", "-n", str(RANKS), "femem", "run", "--estimate", YAML.name]
    assert "--ranks" not in cli

    # (c-page) removed / doubled tags raise
    for tag, lang in ((D1_TAG, "python"), (D2_TAG, "python"), (CLI_TAG, "bash")):
        with pytest.raises(LookupError):
            extract_tagged_block(
                text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag, lang=lang
            )
        doubled = (
            text
            + "\n"
            + f"```{lang} {tag}\n"
            + extract_tagged_block(text, tag, lang=lang)
            + "```\n"
        )
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    # (c-page) d1 naming the tuned file fails (X)
    wrong = d1.replace(f'"{YAML.name}"', '"f_human_tuned_copper.yaml"')
    reason = _x_defect(wrong, D1_TAG)
    _p0(f"[WF-17 s20a] (c-page) d1 naming f_human_tuned_copper.yaml -> {reason}")
    assert reason is not None and "f_human_tuned_copper.yaml" in reason

    # no reading, no artifact-shaped token on the page
    for token in (
        ".json",
        ".npz",
        ".xdmf",
        ".s4p",
        "Study(",
        str(E_RECORD),
        str(N_DOFS_D2_RECORD),
        str(HOLE_CELL_RECORD),
        "GiB",
        "--ranks 2",
    ):
        assert token not in text, token

    # (In) with the build patched to raise; (c-solve) the block resolves the patches
    _scratch(monkeypatch)
    solve_obj, boom = _patch_all_solves(monkeypatch)
    reached = []

    def _raise(*a, **k):
        reached.append(1)
        raise _BuildReached()

    monkeypatch.setattr(FHumanBirdcageGeometry, "build", _raise)
    from fem_em_solver import load_study

    ref = load_study(YAML)
    left_at_1 = d2.replace(D2_LINE, D2_LEFT_AT_1)
    assert left_at_1 != d2
    for tag, block, want_deg in ((D1_TAG, d1, 1), (D2_TAG, d2, 2), ("c-page d2@1", left_at_1, 2)):
        with pytest.raises(_BuildReached):
            _exec_capture(block, tag)
        study = _LAST_NS["study"]
        ok = (
            isinstance(study.geometry, FHumanBirdcageGeometry)
            and study.geometry == ref.geometry
            and study.material == ref.material
            and study.excitation == ref.excitation
            and study.conductor == ref.conductor
            and study.frequency_hz.hex() == ref.frequency_hz.hex()
        )
        _p0(
            f"[WF-17 s20a] (In) {tag}: degree {study.degree} (want {want_deg}); inputs equal "
            f"the tracked file's {ok}; result None {study.result is None}; build reached "
            f"{len(reached)}; resolves patched solves {_resolves_patched(study, solve_obj, boom)}"
        )
        assert ok
        assert study.result is None
        assert _resolves_patched(study, solve_obj, boom)
        if tag == "c-page d2@1":
            assert study.degree != want_deg  # (FH2) on degree fails, by construction
        else:
            assert study.degree == want_deg
    assert len(reached) == 3


_LAST_NS: dict = {}


def _exec_capture(block: str, tag: str) -> dict:
    """``_exec`` that keeps the namespace even when the block raises."""
    _LAST_NS.clear()
    ns = _LAST_NS
    ns["__name__"] = "__main__"
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{tag}", "exec")
    try:
        with contextlib.redirect_stdout(buf):
            exec(code, ns)
    finally:
        _p0(f"[WF-17 s20a] block {tag} stdout (rank 0): {buf.getvalue()!r}")
    return ns


def test_wf17_s20a_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s20a] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- step 26 (no mesh)


def test_wf17_s26_clinical_block_reads_fc2(monkeypatch, capsys, no_build):  # noqa: F811
    from fem_em_solver.cli import main

    comm = MPI.COMM_WORLD
    text = _text()
    block = extract_tagged_block(text, FC_TAG)
    assert f"(../../tests/package/wf13/{FC.name})" in text
    assert sum(ln == "  degree: 1" for ln in FC.read_text().splitlines()) == 1
    assert block.count(D2_LINE) == 1
    for name in (
        "cells",
        "unknowns",
        "degree",
        "ranks",
        "predicted_seconds",
        "predicted_peak_gib",
        "rows_used",
        "not_priced",
    ):
        assert f"est.{name}" in block, name
    assert not _names_used(block) & SOLVE_NAMES

    scratch = SCRATCH / "s26"
    if comm.rank == 0:
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        (scratch / FC.name).write_bytes(FC.read_bytes())
    comm.Barrier()
    monkeypatch.chdir(scratch)

    def _run(src, tag):
        ns = {"__name__": "__main__"}
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            exec(compile(src, f"{PAGE_ID}#{tag}", "exec"), ns)
        return ns, buf.getvalue().strip()

    ns, printed = _run(block, FC_TAG)
    est, study = ns["est"], ns["study"]
    row = FC2_LITERALS
    want = (
        f"{row['cells']} {row['unknowns']} {row['degree']} {row['ranks']} {row['seconds']} "
        f"{row['peak_gib']} ('FC2',) {NOT_PRICED!r}"
    )

    # (Ref) each refusal the page names, on the page's own copy
    with pytest.raises(ValueError, match="xl tier only") as exc:
        study.solve(comm)
    msg = str(exc.value)
    out_dir = scratch / "out"
    capsys.readouterr()
    rc_run = main(["run", FC_D2_NAME, "--output-dir", str(out_dir)])
    err_run = capsys.readouterr().err
    rc_est = main(["run", "--estimate", FC_D2_NAME, "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    names = set(_field_names())
    diff = comm.bcast(_diff(_parse(out, names), est) if comm.rank == 0 else None, root=0)
    rcs_run, rcs_est = comm.allgather(rc_run), comm.allgather(rc_est)

    # Control: the block left at degree 1 reads FC1, not FC2 (runs last: it overwrites the copy)
    left_at_1 = block.replace(D2_LINE, D2_LEFT_AT_1)
    assert left_at_1 != block
    ctrl = _run(left_at_1, "c-page d2@1")[0]["est"]
    ratio = est.unknowns / ctrl.unknowns

    _p0(f"[WF-17 s26] block printed {printed!r}")
    _p0(f"[WF-17 s26] (FC2p) wanted  {want!r}; equal {printed == want}")
    _p0(f"[WF-17 s26] (Ref) solve: {msg}")
    _p0(f"[WF-17 s26] (Ref) femem run: rc {rcs_run}; stderr {err_run.strip()!r}")
    _p0(f"[WF-17 s26] (CLI = API) rc {rcs_est}; differing {diff}")
    _p0(
        f"[WF-17 s26] (c-d1) rows {ctrl.rows_used}, unknowns {ctrl.unknowns}; FC2 / it = "
        f"{ratio:.4f} (UNKNOWN_MISS {UNKNOWN_MISS})"
    )

    assert printed == want
    assert est.rows_used == ("FC2",) and est.not_priced == NOT_PRICED
    assert study.result is None
    assert "<= 3 root-grid nodes" in msg and "xxl tier" in msg and "column set" in msg
    assert rcs_run == [2] * comm.size
    if comm.rank == 0:
        assert "xl tier only" in err_run and "xxl tier" in err_run
    assert rcs_est == [0] * comm.size
    assert diff == set()
    assert not out_dir.exists()
    assert ctrl.rows_used == ("FC1",)
    assert round(ratio, 2) == UNKNOWN_MISS
    assert no_build == []


# ------------------------------------------------------------- build window (-n 2)


def test_wf17_s20a_estimate_blocks_reproduce_fh1_fh2(monkeypatch, capsys):
    if not _knob():
        pytest.skip(
            f"{KNOB} unset: the human-scale estimate window is environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the F-human port Study needs the complex build"
    assert comm.size == RANKS, comm.size
    from fem_em_solver.cli import main

    t0 = time.perf_counter()
    raw = PAGE.read_bytes()
    text = raw.decode()
    d1 = extract_tagged_block(text, D1_TAG)
    d2 = extract_tagged_block(text, D2_TAG)
    _p0(
        f"[WF-17 s20a] page sha256 {_sha(raw)}; {D1_TAG} sha256 {_sha(d1.encode())}; "
        f"{D2_TAG} sha256 {_sha(d2.encode())}"
    )
    _scratch(monkeypatch)
    solve_obj, boom = _patch_all_solves(monkeypatch)

    orig = FHumanBirdcageGeometry.build
    memo, count = {}, {"builds": 0, "hits": 0}

    def _memo_build(self, c=MPI.COMM_WORLD):
        key = (repr(self), c.size)
        if key in memo:
            count["hits"] += 1
            return memo[key]
        count["builds"] += 1
        tb = time.perf_counter()
        memo[key] = orig(self, c)
        _p0(f"[WF-17 s20a] real F-human build {time.perf_counter() - tb:.1f} s")
        return memo[key]

    monkeypatch.setattr(FHumanBirdcageGeometry, "build", _memo_build)
    names = set(_field_names())

    results = {}
    for tag, block, fname in ((D1_TAG, d1, YAML.name), (D2_TAG, d2, D2_NAME)):
        tb = time.perf_counter()
        ns = _exec(block, tag)
        t_block = time.perf_counter() - tb
        est, study = ns["est"], ns["study"]
        out_dir = SCRATCH / f"out-{tag}"
        capsys.readouterr()
        tc = time.perf_counter()
        rc = main(["run", "--estimate", fname, "--output-dir", str(out_dir)])
        t_cli = time.perf_counter() - tc
        out = capsys.readouterr().out
        verdict = None
        if comm.rank == 0:
            parsed = _parse(out, names)
            complete = set(parsed) == names and all(len(v) == 1 for v in parsed.values())
            verdict = dict(complete=complete, differ=_diff(parsed, est) if complete else None)
        else:
            assert _parse(out, names) == {}, out
        verdict = comm.bcast(verdict, root=0)
        rcs = comm.allgather(rc)
        out_exists = comm.allgather(out_dir.exists())
        results[tag] = dict(est=est, study=study, rcs=rcs, verdict=verdict, out_exists=out_exists)
        _p0(f"[WF-17 s20a] {tag}: API estimate {est!r}")
        _p0(
            f"[WF-17 s20a] {tag}: predicted_seconds {est.predicted_seconds:.3f} "
            f"predicted_peak_gib {est.predicted_peak_gib:.4f} at -n {est.ranks} "
            f"(printed only; no -n 2 F-human solve to compare)"
        )
        _p0(
            f"[WF-17 s20a] {tag}: (CLI = API) {fname}: rc {rcs}; complete "
            f"{verdict['complete']}; differing {verdict['differ']}; output dir exists "
            f"{out_exists}; block {t_block:.1f} s, CLI {t_cli:.1f} s"
        )
        _p0(
            f"[WF-17 s20a] {tag}: (NoSolve) result None {study.result is None}; resolves "
            f"patched solves {_resolves_patched(study, solve_obj, boom)}"
        )

    # (c-page) the d2 block left at degree 1 (runs last: it overwrites the copy)
    left_at_1 = d2.replace(D2_LINE, D2_LEFT_AT_1)
    assert left_at_1 != d2
    ctrl = _exec(left_at_1, "c-page d2@1")["est"]
    ctrl_fails = (ctrl.degree != 2, ctrl.unknowns != N_DOFS_D2_RECORD, "FH2" not in ctrl.rows_used)
    _p0(
        f"[WF-17 s20a] (c-page) d2 left at degree 1: degree {ctrl.degree}, unknowns "
        f"{ctrl.unknowns}, rows {ctrl.rows_used} -> (FH2) fails on degree / unknowns / row "
        f"{ctrl_fails}"
    )

    builds = comm.allgather(count["builds"])
    hits = comm.allgather(count["hits"])
    e1, e2 = results[D1_TAG]["est"], results[D2_TAG]["est"]
    _p0(
        f"[WF-17 s20a] (FH1) cells {e1.cells} (HOLE_CELL_RECORD {HOLE_CELL_RECORD}) unknowns "
        f"{e1.unknowns} (E_RECORD {E_RECORD}) degree {e1.degree} rows {e1.rows_used} "
        f"not_priced {e1.not_priced}"
    )
    _p0(
        f"[WF-17 s20a] (FH2) cells {e2.cells} (HOLE_CELL_RECORD {HOLE_CELL_RECORD}) unknowns "
        f"{e2.unknowns} (N_DOFS_D2_RECORD {N_DOFS_D2_RECORD}) degree {e2.degree} rows "
        f"{e2.rows_used} not_priced {e2.not_priced}"
    )
    _p0(
        f"[WF-17 s20a] (B) real builds per rank {builds}, memo hits {hits}; "
        f"node {time.perf_counter() - t0:.1f} s"
    )

    # (FH1)
    assert e1.unknowns == E_RECORD
    assert e1.cells == HOLE_CELL_RECORD
    assert e1.degree == 1
    assert "FH1" in e1.rows_used
    assert e1.not_priced == NOT_PRICED
    # (FH2)
    assert e2.unknowns == N_DOFS_D2_RECORD
    assert e2.cells == HOLE_CELL_RECORD
    assert e2.degree == 2
    assert "FH2" in e2.rows_used
    assert e2.not_priced == NOT_PRICED
    for tag in (D1_TAG, D2_TAG):
        r = results[tag]
        # (NoSolve) / (c-solve)
        assert r["study"].result is None
        assert _resolves_patched(r["study"], solve_obj, boom)
        # (CLI = API)
        assert r["rcs"] == [0] * comm.size
        assert r["verdict"]["complete"]
        assert r["verdict"]["differ"] == set()
        assert r["out_exists"] == [False] * comm.size
    # (c-page) asserted
    assert all(ctrl_fails)
    # (B)
    assert builds == [1] * comm.size
