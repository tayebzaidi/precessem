"""`WF-17` step 20b (§9 item 671, §10 U90b): the loaded-tune page.

``docs/guide/tuning_a_saved_result.md`` carries one tagged block
(``loaded-tune``) that solves the degree-1 PEC hole F-small four-port Study
once, writes its Touchstone file, builds a **fresh** Study with the same
inputs, calls ``load_result`` and ``tune(1)`` and prints ``C_mode``; a second,
**untagged** block shows the F-human command, cited to `WF-16` step 12's
record and **never executed** here.

Brief: docs/planning/chunks/WF-17.md, "§9 item 671".

Anchors: (X) the block's ``C_mode`` vs U94's (CM) record (``U94_CM_RECORD_F``
below) at ``REDUCTION_BAND``, asserted; its ``.9e`` string vs
``ITEM429_C_MODE_F`` printed.  (S0) ``Study.solve`` and ``Study._solve_ports``
are wrapped at class level for the exec and record the Study each call ran on:
exactly one solve, on ``solved``; none on ``fresh``; ``fresh.result`` is a
``LoadedPortResult``.  (D) docrefs 0 dead; only the tagged blocks are executed
(the page holds exactly three python fences, two tagged, the swept one by step
27, §9 item 739; the untagged one is the only one naming
``FHumanBirdcageGeometry`` / ``n_slots`` and is never compiled).  (Ref) ``load_result`` on the F-human Study with a missing file
raises ``FileNotFoundError`` and leaves ``result`` unset (amended by §9 item
685: since `WF-16` step 15 the type check accepts the F-human Study).  (E) the bands block = the
imported constant at ``.0e``.  Negative control (asserted; backed by the same
fixture's C_2 9.873468900e-11 at ``20261015T201110Z_WF-16-s6b-W1.log.gz:940``,
a 0.49 relative distance from C_1): the block with ``fresh.tune(2)`` in place
of ``fresh.tune(1)`` misses (X) by more than ``REDUCTION_BAND``.  Every band
is imported; none is set here.

The solve node runs only with ``FEM_EM_WF17_S20B=1`` (then it asserts the
complex build and ``-n 2``, never skips).  Blocks run in a scratch directory
under the gitignored ``logs/``.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import os
import shutil
import time
from pathlib import Path

import pytest
import yaml
from mpi4py import MPI

from fem_em_solver import Study, StudyConfigError, cli, run_study_file
from fem_em_solver import study as study_mod
from tests.complex_mode import IS_COMPLEX, complex_only
from tests.package.test_wf1_study_yaml import _fresh_dir, _rc_all
from tests.package.test_wf1_yaml_load_result import _loaded_mapping, _no_mesh
from tests.package.test_wf1_yaml_load_result_sweep import ITEM722_C_MODE_F
from tests.package.test_wf16_load_result_sweep import FREQS, SWEEP
from tests.package.test_wf16_study_tune_capacitors import ITEM429_C_MODE_F
from tests.package.test_wf16_study_tune_f_human_copper import _study as _f_human_study
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf17_tune_tutorial import _assert_gate_inputs
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "tuning_a_saved_result.md"
PAGE_ID = "docs/guide/tuning_a_saved_result.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s20b"
S20B_ENV = "FEM_EM_WF17_S20B"
MIN_GUIDE_PAGES = 66

TAG = "loaded-tune"
BANDS_TAG = "loaded-tune-bands"
K1_CALL = "fresh.tune(1)"
K2_CALL = "fresh.tune(2)"

#: U94's (CM) record: the loaded Study's C_mode on this fixture at -n 2,
#: `20261027T111952Z_WF-16-s14-W1.log:950` (printed at .12e, 13 significant digits).
U94_CM_RECORD_F = 6.626796420593e-11
#: `WF-16` step 12's F-human record, cited on the page, never recomputed here:
#: `20261015T213534Z_WF-16-s12-W1.log.gz:7375`, quoted at
#: `tests/package/test_wf1_yaml_tune_f_human.py:75-76`.
F_HUMAN_RECORD_TEXT = "1.181185799070e-10"

# `WF-17` step 27 (§9 item 737): the several-frequency section.
OUT_S27 = REPO / "logs" / "wf17" / "s27"  # gitignored; shared across ranks, not tmp_path
SWEPT_TAG = "swept-tune"
SWEPT_FILE_TAG = "swept-tune-file"
SWEPT_LINES_TAG = "swept-tune-lines"
SWEPT_BANDS_TAG = "swept-tune-bands"
SWEPT_FILE = "f_small_sweep.s4p"
SWEPT_TUPLE = "(1.0e7, 6.4e7, 1.28e8)"
AT_HZ_LINE = "  at_hz: 64000000.0\n"
NEXT_FREQS = (FREQS[1], FREQS[2], FREQS[0])
OUT_OF_SPAN_HZ = 2.0 * FREQS[2]
#: Item 730's (I=) output, `20261101T210812Z_WF-1-s27-W1.log:53,56,59`, in ``FREQS`` order.
ITEM730_HEX = ("0x1.53d78952a75e9p-29", "0x1.0f2f6316ff1dep-34", "0x1.0a3c9bbb69f6ep-36")
SWEPT_REFUSAL_SENTENCES = (
    "An `at_hz` outside the file's span is refused: `tune` raises `ValueError`, and "
    "`femem run` refuses the field `tune.at_hz`, exits with status 2 and writes no output "
    "directory.",
    "`tune` with no `at_hz` raises `ValueError` naming `at_hz`.",
    "`femem run` with no `tune.at_hz` on a swept file refuses that field with status 2.",
)


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _swept_line(f: float, c: float) -> str:
    return f"{f:.6e} Hz: C_mode = {c:.6e} F ({float(c).hex()})"


def _text() -> str:
    return PAGE.read_text()


def _python_fences(text: str) -> list:
    """(opener line, body) of every ```python fence on the page."""
    lines = text.splitlines()
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith("```python"):
            opener, body = lines[i], []
            i += 1
            while lines[i] != "```":
                body.append(lines[i])
                i += 1
            out.append((opener, "\n".join(body) + "\n"))
        i += 1
    return out


def _names(block: str) -> set:
    out = set()
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Name):
            out.add(node.id)
        elif isinstance(node, ast.Attribute):
            out.add(node.attr)
        elif isinstance(node, ast.keyword) and node.arg:
            out.add(node.arg)
        elif isinstance(node, ast.alias):
            out.add(node.name)
    return out


def _bands_line() -> str:
    return f"Saved-file C_mode vs the solved-Study record: {REDUCTION_BAND:.0e}"


def _knob() -> bool:
    raw = os.environ.get(S20B_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S20B_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _exec_block(monkeypatch, block: str, label: str, subdir: str):
    """Rank 0 clears ``SCRATCH_ROOT/subdir``; every rank chdirs there and execs."""
    comm = MPI.COMM_WORLD
    scratch = SCRATCH_ROOT / subdir
    if comm.rank == 0:
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
    comm.Barrier()
    _p0(f"[WF-17 s20b] block {label} sha256 {_sha(block.encode())}")
    monkeypatch.chdir(scratch)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{label}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s20b] block {label} stdout (rank 0): {out!r}")
    return namespace, out


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s20b_blocks_and_refusals(tmp_path):
    text = _text()
    block = extract_tagged_block(text, TAG)
    extract_tagged_block(text, BANDS_TAG, lang="text")
    fences = _python_fences(text)
    openers = [o for o, _ in fences]
    tagged = [b for o, b in fences if o == f"```python {TAG}"]
    untagged = [b for o, b in fences if o == "```python"]
    _p0(f"[WF-17 s20b] (D) python fences {openers}")
    census = [f"```python {TAG}", f"```python {SWEPT_TAG}", "```python"]
    assert openers == census, openers
    # control (§9 item 739): the page with the swept block removed fails the census
    swept = next(b for o, b in fences if o == census[1])
    copy = tmp_path / "page_without_swept_block.md"
    copy.write_text(text.replace(f"{census[1]}\n{swept}```\n", ""), encoding="utf-8")
    control = [o for o, _ in _python_fences(copy.read_text(encoding="utf-8"))]
    _p0(f"[WF-17 s20b] (D) control, swept block removed: python fences {control}")
    assert control == [census[0], census[2]] and control != census, control
    assert tagged == [block]
    fh = untagged[0]
    nb, nf = _names(block), _names(fh)
    _p0(
        f"[WF-17 s20b] tagged block names include load_result {('load_result' in nb)}, tune "
        f"{('tune' in nb)}, solve {('solve' in nb)}, FHuman {('FHumanBirdcageGeometry' in nb)}, "
        f"n_slots {('n_slots' in nb)}; untagged: FHuman {('FHumanBirdcageGeometry' in nf)}, "
        f"n_slots {('n_slots' in nf)}, load_result {('load_result' in nf)}"
    )
    # the tagged block: one solve on `solved`, one load_result + one tune(1) on `fresh`
    assert block.count("Study(") == 1 and block.count(" = f_small_study()") == 2
    assert block.count(".solve(") == 1 and "solved.solve()" in block
    assert block.count(".load_result(") == 1 and "fresh.load_result(path)" in block
    assert block.count(".tune(") == 1 and K1_CALL in block
    assert "write_touchstone(" in block
    assert "FHumanBirdcageGeometry" not in nb and "n_slots" not in nb
    # the untagged block is the F-human command and never loads
    assert {"FHumanBirdcageGeometry", "n_slots", "LeontovichWall"} <= nf
    assert "load_result" not in nf and "tune(1, n_slots=16)" in fh
    assert F_HUMAN_RECORD_TEXT in text
    assert "20261015T213534Z_WF-16-s12-W1.log.gz" in text
    assert "tests/package/test_wf1_yaml_tune_f_human.py" in text
    assert "(tuning_a_study.md)" in text
    for token in ("AED", "HFSS"):
        assert token not in text, token

    # (c-page), by construction: tag removed / doubled raise
    for tag, lang in ((TAG, "python"), (BANDS_TAG, "text")):
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
    _p0("[WF-17 s20b] (c-page) tags removed / doubled raise LookupError for both blocks")

    # (Ref) amended by §9 item 685 (`WF-16` step 15): load_result accepts the F-human
    # Study now, so a missing file is refused by the file check, not the type check
    missing = tmp_path / "no_such_file.s32p"
    st = _f_human_study()
    with pytest.raises(FileNotFoundError) as ei:
        st.load_result(missing)
    assert st.result is None
    _p0(f"[WF-17 s20b] (Ref) F-human load_result on a missing file: FileNotFoundError {ei.value}")


def test_wf17_s20b_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    got = block.strip("\n").splitlines()
    _p0(f"[WF-17 s20b] (E) bands block {got} (rendered {[_bands_line()]})")
    assert got == [_bands_line()]
    altered = block.replace(f"{REDUCTION_BAND:.0e}", f"{REDUCTION_BAND * 10:.0e}", 1)
    assert altered.strip("\n").splitlines() != [_bands_line()]


def test_wf17_s20b_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s20b] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s20b_loaded_tune_block(monkeypatch):
    if not _knob():
        pytest.skip(
            f"{S20B_ENV} unset: the loaded-tune solve window is environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port Study needs the complex build"
    assert comm.size == 2, comm.size
    block = extract_tagged_block(_text(), TAG)
    assert block.count(K1_CALL) == 1
    k2_block = block.replace(K1_CALL, K2_CALL, 1)

    real_solve, real_ports = Study.solve, Study._solve_ports

    def _run(src, label):
        calls = []

        def solve_spy(self, *a, **k):
            calls.append(("solve", self))
            return real_solve(self, *a, **k)

        def ports_spy(self, *a, **k):
            calls.append(("_solve_ports", self))
            return real_ports(self, *a, **k)

        t0 = time.perf_counter()
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(Study, "solve", solve_spy)
            mp.setattr(Study, "_solve_ports", ports_spy)
            ns, out = _exec_block(monkeypatch, src, label, label)
        return ns, out, calls, time.perf_counter() - t0

    ns, out, calls, t1 = _run(block, "k1")
    ns2, out2, calls2, t2 = _run(k2_block, "k2")

    # (S0)
    solved, fresh = ns["solved"], ns["fresh"]
    on_solved = [n for n, s in calls if s is solved]
    on_fresh = [n for n, s in calls if s is fresh]
    _p0(
        f"[WF-17 s20b] (S0) calls {[n for n, _ in calls]}; on solved {on_solved}; on fresh "
        f"{on_fresh}; fresh.result {type(fresh.result).__name__}; block {t1:.1f} s"
    )
    # (In) both Studies are the gate's inputs
    _assert_gate_inputs(solved, "solved")
    _assert_gate_inputs(fresh, "fresh")

    # (X)
    tuned = ns["tuned"]
    c = float(tuned.c_mode_farad)
    rel = abs(c / U94_CM_RECORD_F - 1.0)
    rel_g = comm.allreduce(rel, op=MPI.MAX)
    _p0(
        f"[WF-17 s20b] (X) block k {tuned.k}, C_mode {c:.12e} F (hex {c.hex()}) vs U94 (CM) "
        f"{U94_CM_RECORD_F:.12e}: rel {rel_g:.3e} (ASSERTED <= REDUCTION_BAND "
        f"{REDUCTION_BAND:.0e}); '{c:.9e}' vs ITEM429_C_MODE_F '{ITEM429_C_MODE_F:.9e}' "
        f"equal {f'{c:.9e}' == f'{ITEM429_C_MODE_F:.9e}'} (printed)"
    )

    # negative control: the block with tune(2)
    c2 = float(ns2["tuned"].c_mode_farad)
    miss = abs(c2 / U94_CM_RECORD_F - 1.0)
    miss_g = comm.allreduce(miss, op=MPI.MIN)
    on_fresh2 = [n for n, s in calls2 if s is ns2["fresh"]]
    _p0(
        f"[WF-17 s20b] (c-k2) tune(2) block k {ns2['tuned'].k}, C {c2:.12e} F: miss {miss_g:.4e} "
        f"({miss_g / REDUCTION_BAND:.1f} x REDUCTION_BAND; ASSERTED > {REDUCTION_BAND:.0e}); "
        f"solves on its fresh Study {on_fresh2}; block {t2:.1f} s"
    )

    assert on_solved == ["solve", "_solve_ports"], on_solved
    assert on_fresh == []
    assert len(calls) == 2
    assert isinstance(fresh.result, study_mod.LoadedPortResult)
    assert tuned.k == 1
    assert rel_g <= REDUCTION_BAND
    assert ns2["tuned"].k == 2 and on_fresh2 == []
    assert miss_g > REDUCTION_BAND


# ------------------------------------------------------------- `WF-17` step 27 (§9 item 737)
#
# The page's "A file holding several frequencies" section, on U71's file (gitignored;
# absent ⇒ fail, never re-solve).  (Lines) the block's stdout and the study file's
# summaries at each ``at_hz`` hold the page's three lines; (C=) each value against item
# 722's record column at ``REDUCTION_BAND``; (Sel) the selected frequency exactly; (Hex)
# the value's ``float.hex`` equal to item 730's (I=) output; (S0) 0 solves.  Control: both
# routes with every ``at_hz`` moved one frequency along miss by > ``REDUCTION_BAND``.
# (Ref) the page's refusal sentence, on the page's own block and study file.


@complex_only
def test_wf17_s27_swept_file_section(capsys, monkeypatch):
    comm = MPI.COMM_WORLD
    assert SWEEP.is_file(), f"{SWEEP} is absent: this test never re-solves"
    text = _text()
    block = extract_tagged_block(text, SWEPT_TAG)
    file_text = extract_tagged_block(text, SWEPT_FILE_TAG, lang="yaml")
    lines = extract_tagged_block(text, SWEPT_LINES_TAG, lang="text").strip("\n").splitlines()
    bands = extract_tagged_block(text, SWEPT_BANDS_TAG, lang="text").strip("\n").splitlines()
    assert block.count(SWEPT_TUPLE) == 1 and file_text.count(AT_HZ_LINE) == 1
    root = OUT_S27 / f"run_n{comm.size}"
    _fresh_dir(root, comm)
    if comm.rank == 0:
        root.mkdir(parents=True)
        shutil.copy(SWEEP, root / SWEPT_FILE)
    comm.Barrier()
    monkeypatch.chdir(root)

    # (E), (In): the band is the imported constant; the file is the tracked study + two keys
    expected_file = _loaded_mapping(SWEPT_FILE)
    expected_file["tune"] = {"k": 1, "at_hz": FREQS[1]}
    assert yaml.safe_load(file_text) == expected_file
    band_ok = bands == [
        f"Swept-file C_mode vs the single-frequency tune of its column: " f"{REDUCTION_BAND:.0e}"
    ]

    def exec_block(src, label):
        buf = io.StringIO()
        code = compile(src, f"{PAGE_ID}#{label}", "exec")
        ns = {"__name__": "__main__"}
        with contextlib.redirect_stdout(buf):
            exec(code, ns)
        return ns, buf.getvalue()

    def run_file(src, label, f):
        yml = root / f"{label}.yaml"
        if comm.rank == 0:
            yml.write_text(src.replace(AT_HZ_LINE, f"  at_hz: {f!r}\n"))
        comm.Barrier()
        run = run_study_file(yml, output_dir=root / label, comm=comm)
        capsys.readouterr()
        rc = cli.main(["run", str(yml), "--output-dir", str(root / f"{label}_cli")])
        capsys.readouterr()
        return run.summary["tune"], rc

    def refused(src, label, f):
        yml = root / f"{label}.yaml"
        if comm.rank == 0:
            yml.write_text(src.replace(AT_HZ_LINE, f"  at_hz: {f!r}\n"))
        comm.Barrier()
        out = root / label
        with pytest.raises(StudyConfigError) as err:
            run_study_file(yml, output_dir=out, comm=comm)
        capsys.readouterr()
        rc = cli.main(["run", str(yml), "--output-dir", str(out)])
        capsys.readouterr()
        return err.value.field, _rc_all(rc, comm), out.exists()

    next_tuple = "(" + ", ".join(repr(f) for f in NEXT_FREQS) + ")"
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        ns, out = exec_block(block, "swept")
        ns_next, out_next = exec_block(block.replace(SWEPT_TUPLE, next_tuple), "swept-next")
        file_runs = [run_file(file_text, f"f{i}", f) for i, f in enumerate(FREQS)]
        next_runs = [run_file(file_text, f"n{i}", f) for i, f in enumerate(NEXT_FREQS)]
        # (Ref) the page's sentence, on the page's own copy
        with pytest.raises(ValueError, match="outside the file's span") as e_span:
            exec_block(block.replace(SWEPT_TUPLE, f"({OUT_OF_SPAN_HZ!r},)"), "span")
        with pytest.raises(ValueError, match="at_hz") as e_none:
            exec_block(block.replace("study.tune(1, at_hz=f)", "study.tune(1)"), "none")
        ref_file = refused(file_text, "ref", OUT_OF_SPAN_HZ)
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)

    # Printed before any assert, so a red window still carries every reading.
    block_lines = out.splitlines()
    rows = [(float(t.frequency_hz), float(t.c_mode_farad)) for t in ns["rows"]]
    summ = [(float(b["frequency_hz"]), float(b["c_mode_farad"])) for b, _rc in file_runs]
    summ_lines = [_swept_line(f, c) for f, c in summ]
    rcs = _rc_all(max(abs(rc) for _b, rc in file_runs), comm)
    _p0(
        f"\n[WF-17 s27] block sha256 {_sha(block.encode())}, file sha256 "
        f"{_sha(file_text.encode())}; (S0) mesh and solve patched: {n_calls} calls"
    )
    bad = []
    for i, f in enumerate(FREQS):
        sel, c = rows[i]
        rel = abs(c / ITEM722_C_MODE_F[i] - 1.0)
        rel_f = abs(summ[i][1] / ITEM722_C_MODE_F[i] - 1.0)
        gate = (
            lines[i] in block_lines,
            summ_lines[i] == lines[i],
            sel == f,
            summ[i][0] == f,
            c.hex() == ITEM730_HEX[i],
            summ[i][1].hex() == ITEM730_HEX[i],
            rel <= REDUCTION_BAND,
            rel_f <= REDUCTION_BAND,
        )
        _p0(
            f"[WF-17 s27 lines] {f:.1e} Hz: page {lines[i]!r}; in block stdout {gate[0]}, equals "
            f"summary line {gate[1]}"
        )
        _p0(
            f"[WF-17 s27 C=] {f:.1e} Hz: selected block {sel:.1e} {gate[2]}, file {summ[i][0]:.1e} "
            f"{gate[3]}; hex block {gate[4]}, file {gate[5]}; rel to item 722 {rel:.3e} / "
            f"{rel_f:.3e} (REDUCTION_BAND {REDUCTION_BAND:.0e})"
        )
        if not all(gate):
            bad.append(f)
    rows_n = [float(t.c_mode_farad) for t in ns_next["rows"]]
    next_n = [float(b["c_mode_farad"]) for b, _rc in next_runs]
    miss = min(abs(c / rec - 1.0) for c, rec in zip(rows_n, ITEM722_C_MODE_F))
    miss_f = min(abs(c / rec - 1.0) for c, rec in zip(next_n, ITEM722_C_MODE_F))
    _p0(
        f"[WF-17 s27 ctl] every at_hz one frequency along: least miss block {miss:.6f}, file "
        f"{miss_f:.6f} (asserted > {REDUCTION_BAND:.0e}); stdout differs "
        f"{out_next != out}"
    )
    _p0(f"[WF-17 s27 Ref] block, out of span: {e_span.value}")
    _p0(f"[WF-17 s27 Ref] block, no at_hz: {e_none.value}")
    _p0(f"[WF-17 s27 Ref] file, out of span: field, rc, output dir {ref_file}; run rc {rcs}")

    # The sentence the (Ref) above backs is the page's, whitespace normalised.
    prose = " ".join(text.split())
    for sentence in SWEPT_REFUSAL_SENTENCES:
        assert sentence in prose, sentence
    assert "AED" not in text and "HFSS" not in text
    assert n_calls == 0
    assert band_ok, bands
    assert bad == []
    assert rcs == (0, 0)
    assert miss > REDUCTION_BAND and miss_f > REDUCTION_BAND
    assert ref_file == ("tune.at_hz", (2, 2), False)
