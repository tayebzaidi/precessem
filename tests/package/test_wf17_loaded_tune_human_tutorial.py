"""`WF-17` step 24 (§9 item 691, §10 U103): the saved human-scale tune page.

``docs/guide/tuning_a_saved_human_scale_result.md`` carries one tagged block
(``loaded-tune-human``) that solves the copper F-human 32-port Study once,
writes its Touchstone file, builds a **fresh** Study with the same inputs,
calls ``load_result`` and ``tune(1, n_slots=16)`` and prints ``C_mode``.

Brief: docs/planning/chunks/WF-17.md, "§9 item 691".

Anchors: (X) the block's ``C_mode`` vs ``WF12_STEP3_C_MODE_CU_F`` at
``REDUCTION_BAND`` (U100's (CM), rel 3.421e-13 at
``20261028T004356Z_WF-16-s15-W1.log:7353``), asserted; the page's printed
``.4e`` equals the record's ``.4e``.  (S0) class-level spies on
``Study.solve`` / ``Study._solve_ports``: one call each, on ``solved``; none
on ``fresh``.  (In) both Studies equal ``_study()`` of
``test_wf16_study_tune_f_human_copper.py`` field by field.  (D) docrefs 0
dead; exactly one python fence on the page, tagged.  (E) the bands block =
the imported constant at ``.0e``.  Negative control (asserted): the page's
written file with ports 1 and 2 swapped (U100's ``_swap12_n`` /
``_write_named``) misses the record by more than ``LEG_D0_REPRODUCTION_BAND``
(U100's record 0.4254907, ``:7356``); a ``ValueError`` from ``tune`` counts
as a miss.  Every band is imported; none is set here.

The solve node runs only with ``FEM_EM_WF17_S24=1`` (unset => skip; any other
value raises); then it asserts the complex build and ``-n 8``.  Blocks run in
a scratch directory under the gitignored ``logs/``.
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
from mpi4py import MPI

from fem_em_solver import Study
from fem_em_solver import study as study_mod
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_load_result_f_human import _swap12_n, _write_named
from tests.package.test_wf16_study_tune_f_human_copper import (
    RECORD_RANK_WIDTH,
    WF12_STEP3_C_MODE_CU_F,
)
from tests.package.test_wf16_study_tune_f_human_copper import _study as _f_human_study
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "tuning_a_saved_human_scale_result.md"
PAGE_ID = "docs/guide/tuning_a_saved_human_scale_result.md"
OLD_PAGE = REPO / "docs" / "guide" / "tuning_a_saved_result.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s24"
S24_ENV = "FEM_EM_WF17_S24"
MIN_GUIDE_PAGES = 66

TAG = "loaded-tune-human"
BANDS_TAG = "loaded-tune-human-bands"
K1_CALL = "fresh.tune(1, n_slots=16)"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _python_openers(text: str) -> list:
    return [ln for ln in text.splitlines() if ln.startswith("```python")]


def _bands_line() -> str:
    return f"Saved-file C_mode vs the solved-Study record: {REDUCTION_BAND:.0e}"


def _knob() -> bool:
    raw = os.environ.get(S24_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S24_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _factory_only(block: str) -> str:
    """The block's imports and its ``def f_human_study`` only (no solve)."""
    tree = ast.parse(block)
    keep = [
        n
        for n in tree.body
        if isinstance(n, (ast.Import, ast.ImportFrom))
        or (isinstance(n, ast.FunctionDef) and n.name == "f_human_study")
    ]
    return ast.unparse(ast.Module(body=keep, type_ignores=[]))


def _assert_inputs(st, label):
    ref = _f_human_study()
    _p0(
        f"[WF-17 s24] (In) {label}: geometry {st.geometry!r}; material {st.material!r}; "
        f"excitation {st.excitation!r}; frequency {float(st.frequency_hz).hex()} (want "
        f"{float(ref.frequency_hz).hex()}); degree {st.degree}; conductor {st.conductor!r}"
    )
    assert type(st.geometry) is type(ref.geometry) and st.geometry == ref.geometry
    assert st.material == ref.material
    assert float(st.material.sigma).hex() == float(ref.material.sigma).hex()
    assert float(st.material.epsilon_r).hex() == float(ref.material.epsilon_r).hex()
    assert type(st.excitation) is type(ref.excitation) and st.excitation == ref.excitation
    assert float(st.frequency_hz).hex() == float(ref.frequency_hz).hex()
    assert st.degree == ref.degree
    assert type(st.conductor) is type(ref.conductor) and st.conductor == ref.conductor


def _exec_block(monkeypatch, block: str, label: str, subdir: str):
    """Rank 0 clears ``SCRATCH_ROOT/subdir``; every rank chdirs there and execs."""
    comm = MPI.COMM_WORLD
    scratch = SCRATCH_ROOT / subdir
    if comm.rank == 0:
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
    comm.Barrier()
    _p0(f"[WF-17 s24] block {label} sha256 {_sha(block.encode())}")
    monkeypatch.chdir(scratch)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{label}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s24] block {label} stdout (rank 0): {out!r}")
    return namespace, out, scratch


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s24_block_shape_and_inputs(monkeypatch):
    for bad in ("2", "yes"):
        monkeypatch.setenv(S24_ENV, bad)
        with pytest.raises(ValueError, match=S24_ENV):
            _knob()
    monkeypatch.delenv(S24_ENV, raising=False)
    assert _knob() is False

    text = _text()
    block = extract_tagged_block(text, TAG)
    extract_tagged_block(text, BANDS_TAG, lang="text")
    openers = _python_openers(text)
    _p0(f"[WF-17 s24] (D) python fences {openers}")
    assert openers == [f"```python {TAG}"], openers
    assert block.count("Study(") == 1 and block.count(" = f_human_study()") == 2
    assert block.count(".solve(") == 1 and "solved.solve()" in block
    assert block.count(".load_result(") == 1 and "fresh.load_result(path)" in block
    assert block.count(".tune(") == 1 and K1_CALL in block
    assert block.count("write_touchstone(") == 1
    for token in ("AED", "HFSS"):
        assert token not in text, token
    # the old page links here, in prose
    old = OLD_PAGE.read_text()
    assert "(tuning_a_saved_human_scale_result.md)" in old
    assert "(tuning_a_saved_result.md)" in text

    # (In), no mesh: the block's factory builds the gate's Study
    ns = {"__name__": "wf17_s24_factory"}
    exec(compile(_factory_only(block), f"{PAGE_ID}#factory", "exec"), ns)
    _assert_inputs(ns["f_human_study"](), "factory (W0)")

    # tag removed / doubled raise
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
    _p0("[WF-17 s24] (c-page) tags removed / doubled raise LookupError for both blocks")


def test_wf17_s24_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    got = block.strip("\n").splitlines()
    _p0(f"[WF-17 s24] (E) bands block {got} (rendered {[_bands_line()]})")
    assert got == [_bands_line()]
    altered = block.replace(f"{REDUCTION_BAND:.0e}", f"{REDUCTION_BAND * 10:.0e}", 1)
    assert altered.strip("\n").splitlines() != [_bands_line()]


def test_wf17_s24_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s24] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 8)


def test_wf17_s24_loaded_tune_human_block(monkeypatch):
    if not _knob():
        pytest.skip(f"{S24_ENV} unset: the copper F-human solve is a heavy -n 8 window")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the 32-port Study needs the complex build"
    assert comm.size == RECORD_RANK_WIDTH, (comm.size, RECORD_RANK_WIDTH)
    block = extract_tagged_block(_text(), TAG)

    real_solve, real_ports = Study.solve, Study._solve_ports
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
        ns, out, scratch = _exec_block(monkeypatch, block, "k1", "k1")
    t1 = time.perf_counter() - t0

    # (S0)
    solved, fresh = ns["solved"], ns["fresh"]
    on_solved = [n for n, s in calls if s is solved]
    on_fresh = [n for n, s in calls if s is fresh]
    _p0(
        f"[WF-17 s24] (S0) calls {[n for n, _ in calls]}; on solved {on_solved}; on fresh "
        f"{on_fresh}; fresh.result {type(fresh.result).__name__}; block {t1:.1f} s"
    )
    # (In)
    _assert_inputs(solved, "solved")
    _assert_inputs(fresh, "fresh")

    # (X)
    tuned = ns["tuned"]
    c = float(tuned.c_mode_farad)
    rel = abs(c / WF12_STEP3_C_MODE_CU_F - 1.0)
    rel_g = comm.allreduce(rel, op=MPI.MAX)
    printed = f"C_mode = {WF12_STEP3_C_MODE_CU_F:.4e} F"
    _p0(
        f"[WF-17 s24] (X) block k {tuned.k}, C_mode {c:.12e} F (hex {c.hex()}) vs "
        f"WF12_STEP3_C_MODE_CU_F {WF12_STEP3_C_MODE_CU_F:.12e}: rel {rel_g:.3e} (ASSERTED <= "
        f"REDUCTION_BAND {REDUCTION_BAND:.0e}); printed '{printed}' in stdout "
        f"{printed in out}"
    )

    # negative control: the page's file with ports 1 <-> 2 swapped (data permuted, names kept)
    loaded = fresh.result
    sw = _write_named(
        scratch / "swap12.s32p",
        fresh.frequency_hz,
        _swap12_n(loaded.s_matrix),
        list(loaded.port_ids),
    )
    st_sw = _f_human_study()
    st_sw.load_result(sw)
    try:
        cw = float(st_sw.tune(1, n_slots=16, comm=comm).c_mode_farad)
        miss = abs(cw / WF12_STEP3_C_MODE_CU_F - 1.0)
        miss_g = comm.allreduce(miss, op=MPI.MIN)
        _p0(
            f"[WF-17 s24] (c-swap) P1<->P2 swapped: C {cw:.12e} F, miss {miss_g:.7f} "
            f"(ASSERTED > LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e})"
        )
    except ValueError as exc:  # a non-capacitive read is a miss too; printed
        miss_g = float("inf")
        _p0(f"[WF-17 s24] (c-swap) P1<->P2 swapped: tune raised {exc} (miss = inf)")

    assert on_solved == ["solve", "_solve_ports"], on_solved
    assert on_fresh == []
    assert len(calls) == 2
    assert isinstance(fresh.result, study_mod.LoadedPortResult)
    assert len(loaded.port_ids) == 32
    assert tuned.k == 1
    assert rel_g <= REDUCTION_BAND
    assert printed in out
    assert miss_g > LEG_D0_REPRODUCTION_BAND
