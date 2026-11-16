"""`WF-17` step 10: the estimate-before-a-solve tutorial page reproduces the
F-small degree-1 `WF-13` estimate through the API block and ``femem run
--estimate``, with no solve.

Brief: docs/planning/chunks/WF-17.md, "§9 item 602".  The mesh node runs only
with ``FEM_EM_WF17_S10=1`` (and then asserts the complex build and ``-n 2``,
never skips).  The other nodes build no mesh.

Anchors: (X) one block per tag, one ``load_study(`` call, one ``.estimate(``
call, no solve name, the file argument read from the AST, the CLI line's tokens;
(F) the fields block = the ``ResourceEstimate`` field names in order; (U) cells
= ``ANS4_EXACT_STRIP_CELL_RECORD``, unknowns = ``PRICE_TABLE`` P3; (P) the
prediction = ``_predict`` bitwise; (K) no solve; (CLI = API) field by field;
(D) 0 dead links.  Controls (rule (e)): (c-deg) asserted, backed; (c-width)
and (c-page) by construction.  Ratios to P3 and build times are printed, never
asserted (P3 is a fitted row).
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

from fem_em_solver.study import PRICE_TABLE, ResourceEstimate, _predict

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf13_cli_estimate import _diff, _parse, _patch_solves
from tests.package.test_wf13_study_estimate import F_SMALL_DEG2_UNKNOWNS_RECORD
from tests.package.test_wf17_first_study_tutorial import _checker_stdout, extract_tagged_block
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "estimate_before_a_solve.md"
PAGE_ID = "docs/guide/estimate_before_a_solve.md"
YAML = REPO / "tests" / "package" / "wf1" / "f_small_four_port.yaml"
SCRATCH = REPO / "logs" / "wf17" / "s10"
S10_ENV = "FEM_EM_WF17_S10"
MIN_GUIDE_PAGES = 59
RANKS = 2

PY_TAG = "estimate-python"
CLI_TAG = "estimate-cli"
FIELDS_TAG = "estimate-fields"

LINK_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _p3() -> dict:
    rows = [r for r in PRICE_TABLE if r["row"] == "P3"]
    assert len(rows) == 1
    return rows[0]


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


def _attr_calls(block: str, name: str) -> list:
    return [
        n for n in ast.walk(ast.parse(block))
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == name
    ]


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


def _file_arg(block: str) -> str:
    calls = _calls_named(block, "load_study")
    if len(calls) != 1:
        raise LookupError(f"{len(calls)} load_study calls; need exactly one")
    return ast.literal_eval(calls[0].args[0])


def _fields_defect(block: str):
    """None when the fields block equals the dataclass field names in order, else the reason."""
    got = [ln.strip() for ln in block.strip("\n").splitlines()]
    want = _field_names()
    if got == want:
        return None
    missing = [w for w in want if w not in got]
    extra = [g for g in got if g not in want]
    return f"fields {got} != {want}; missing {missing}; extra {extra}"


def _knob() -> bool:
    raw = os.environ.get(S10_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S10_ENV}={raw!r}: only unset/empty or '1' are accepted")


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s10_blocks_and_fields(tmp_path):
    text = _text()
    block = extract_tagged_block(text, PY_TAG)
    n_load = len(_calls_named(block, "load_study"))
    n_est = len(_attr_calls(block, "estimate"))
    used = _names_used(block)
    name = _file_arg(block)
    _p0(f"[WF-17 s10] (X) {PY_TAG}: load_study calls {n_load}, .estimate calls {n_est}; "
        f"file {name!r} (want {YAML.name!r}); solve names "
        f"{sorted(used & {'solve', 'run_study_file'})}")
    assert block.count("load_study(") == 1 and n_load == 1
    assert block.count(".estimate(") == 1 and n_est == 1
    assert not (used & {"solve", "run_study_file"})
    assert ".solve(" not in block and "run_study_file" not in block
    assert name == YAML.name
    assert YAML.is_file()
    assert f"(../../tests/package/wf1/{YAML.name})" in text

    cli = extract_tagged_block(text, CLI_TAG, lang="bash")
    lines = [ln for ln in cli.splitlines() if ln.strip()]
    _p0(f"[WF-17 s10] (X) {CLI_TAG}: {lines!r}")
    assert len(lines) == 1
    toks = lines[0].split()
    assert toks[:3] == ["mpiexec", "-n", str(RANKS)]
    assert toks[3:] == ["femem", "run", "--estimate", YAML.name]

    # (F)
    fields = extract_tagged_block(text, FIELDS_TAG, lang="text")
    _p0(f"[WF-17 s10] (F) fields block {fields.split()} (dataclass {_field_names()})")
    assert _fields_defect(fields) is None, _fields_defect(fields)

    # (c-page), by construction: removed / doubled tags raise
    for tag, lang in ((PY_TAG, "python"), (CLI_TAG, "bash"), (FIELDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    _p0("[WF-17 s10] (c-page) removed / doubled tags raise LookupError for all three blocks")

    # (c-page), by construction: a tmp copy with rows_used dropped fails (F), naming it
    altered_page = tmp_path / "estimate_before_a_solve.md"
    altered_page.write_text(text.replace("rows_used\n", "", 1))
    reason = _fields_defect(extract_tagged_block(altered_page.read_text(), FIELDS_TAG, lang="text"))
    _p0(f"[WF-17 s10] (c-page) rows_used dropped -> {reason}")
    assert reason is not None and "rows_used" in reason

    # no reading, no artifact-shaped token on the page
    p3 = _p3()
    for token in (".json", ".npz", ".xdmf", ".s4p", "Study(", str(p3["unknowns"]),
                  str(p3["cells"]), str(ANS4_EXACT_STRIP_CELL_RECORD), "GiB"):
        assert token not in text, token


def test_wf17_s10_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s10] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- mesh window (-n 2)


def test_wf17_s10_estimate_block_reproduces_p3(monkeypatch, capsys):
    if not _knob():
        pytest.skip(f"{S10_ENV} unset: the tutorial estimate window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port port Study needs the complex build"
    assert comm.size == RANKS, comm.size
    from fem_em_solver.cli import main

    t0 = time.perf_counter()
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), PY_TAG)
    name = _file_arg(block)
    if comm.rank == 0:
        print(f"[WF-17 s10] page sha256 {_sha(raw)}; block {PY_TAG} sha256 "
              f"{_sha(block.encode())}", flush=True)
        if SCRATCH.exists():
            shutil.rmtree(SCRATCH)
        SCRATCH.mkdir(parents=True)
        data = YAML.read_bytes()
        (SCRATCH / name).write_bytes(data)
        print(f"[WF-17 s10] copied {name} sha256 {_sha(data)} (source sha256 "
              f"{_sha(YAML.read_bytes())})", flush=True)
    comm.Barrier()
    monkeypatch.chdir(SCRATCH)

    # (K) the solve paths raise if reached
    _patch_solves(monkeypatch)

    ns = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{PY_TAG}", "exec")
    t_b = time.perf_counter()
    with contextlib.redirect_stdout(buf):
        exec(code, ns)
    t_block = time.perf_counter() - t_b
    _p0(f"[WF-17 s10] block {PY_TAG} stdout (rank 0): {buf.getvalue()!r}")
    est = ns["est"]
    study = ns["study"]

    # (CLI = API), in-process, after the block
    out_dir = SCRATCH / "out"
    capsys.readouterr()
    t_c = time.perf_counter()
    rc = main(["run", "--estimate", name, "--output-dir", str(out_dir)])
    t_cli = time.perf_counter() - t_c
    out = capsys.readouterr().out
    names = _field_names()
    verdict = None
    if comm.rank == 0:
        parsed = _parse(out, names)
        complete = set(parsed) == set(names) and all(len(v) == 1 for v in parsed.values())
        differ = _diff(parsed, est) if complete else None
        verdict = dict(complete=complete, differ=differ)
        print("[WF-17 s10] printed by `femem run --estimate`:\n" + out.rstrip(), flush=True)
    else:
        assert _parse(out, names) == {}, out
    verdict = comm.bcast(verdict, root=0)
    rcs = comm.allgather(rc)
    out_exists = comm.allreduce(int(out_dir.exists()), op=MPI.MAX)

    p3 = _p3()
    pred = _predict(est.unknowns, 1, RANKS)
    pred_w1 = _predict(est.unknowns, 1, 1)
    deg2 = _predict(F_SMALL_DEG2_UNKNOWNS_RECORD, 2, RANKS)
    width_ratio = pred_w1[1] / est.predicted_seconds

    _p0(f"[WF-17 s10] API estimate: {est!r}")
    _p0(f"[WF-17 s10] (U) cells {est.cells} (record ANS4_EXACT_STRIP_CELL_RECORD "
        f"{ANS4_EXACT_STRIP_CELL_RECORD}; predicted 117206); unknowns {est.unknowns} (P3 "
        f"{p3['unknowns']}; predicted 140460); degree {est.degree}; ranks {est.ranks}")
    _p0(f"[WF-17 s10] (P) estimate ({est.predicted_peak_gib!r}, {est.predicted_seconds!r}, "
        f"{est.rows_used!r}) vs _predict {pred!r}")
    _p0(f"[WF-17 s10] not_priced {est.not_priced!r} (predicted ())")
    _p0(f"[WF-17 s10] (K) study.result is None {study.result is None}")
    _p0(f"[WF-17 s10] (CLI = API) rc per rank {rcs}; complete {verdict['complete']}; differing "
        f"{verdict['differ']}; {out_dir.name} exists {bool(out_exists)}")
    _p0(f"[WF-17 s10] (c-deg) deg-2 unknowns {F_SMALL_DEG2_UNKNOWNS_RECORD} vs {est.unknowns}; "
        f"_predict deg 2 {deg2[:2]!r} vs {(est.predicted_peak_gib, est.predicted_seconds)!r}")
    _p0(f"[WF-17 s10] (c-width) _predict(N, 1, 1) seconds / estimate seconds = {width_ratio!r}")
    _p0(f"[WF-17 s10] printed, never asserted: predicted_seconds / P3 seconds = "
        f"{est.predicted_seconds / p3['seconds']:.4f}; predicted_peak_gib / P3 peak_gib = "
        f"{est.predicted_peak_gib / p3['peak_gib']:.4f} (in-sample; P3 is a fitted row)")
    _p0(f"[WF-17 s10] build times: block (mesh 1) {t_block:.2f} s; CLI (mesh 2) {t_cli:.2f} s; "
        f"node {time.perf_counter() - t0:.1f} s")

    # (U)
    assert est.cells == ANS4_EXACT_STRIP_CELL_RECORD
    assert est.unknowns == p3["unknowns"]
    assert est.degree == 1
    assert est.ranks == RANKS
    # (P)
    assert (est.predicted_peak_gib, est.predicted_seconds, est.rows_used) == pred
    # (K)
    assert study.result is None
    # (CLI = API)
    assert rcs == [0] * comm.size
    assert verdict["complete"]
    assert verdict["differ"] == set()
    assert not out_exists
    # (c-deg)
    assert F_SMALL_DEG2_UNKNOWNS_RECORD != est.unknowns
    assert deg2[:2] != (est.predicted_peak_gib, est.predicted_seconds)
    # (c-width)
    assert abs(width_ratio - 2.0) <= 1e-15 * 2.0
