"""`WF-17` step 1 (§9 item 257, chain U7): the first tutorial is a test.

The committed page `docs/guide/first_study.md` carries one fenced block tagged
```` ```python first-study ````. This module reads the committed page, extracts
that block verbatim and ``exec``s it (degree 2 by default, ``-n 2``). Asserted
(brief: docs/planning/chunks/WF-17.md#§9-item-257):

(X)   exactly one tagged block, with exactly one ``Study(`` call;
(K)   the block configures `WF-15`'s Study exactly (no solve), with no
      ``degree=`` keyword;
(E)   the page's expected-output block states `WF-15`'s record exactly;
(D)   ``docrefs --docs-root docs/guide --check-links`` reads dead == 0;
(Deg) the solve is degree 2 by default, structurally;
(R)   the block's mean SAR reproduces `WF-15`'s degree-2 record within
      ``L0_BAND``;
(G)   inside `MAT-4`'s imported band against the closed form;
(Out) the block prints what it returns.

Negative controls: (c1) untagged / doubly-tagged page raises, (c2) a changed σ
is caught by (K)'s comparator, (c3) a planted dead link reads dead == 1,
(c4) a planted ``degree=1`` block misses (R) and lands on `MAT-4`'s degree-1
record — all asserted.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest
from mpi4py import MPI

from fem_em_solver import (
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    Study,
)

from tests.complex_mode import complex_only
from tests.materials.test_phantom_material_model import MAT4_BAND
from tests.package.test_wf15_study_lossy_sphere import (
    GEOMETRY,
    MAT4_DEG1_RECORD,
    _closed,
)
from tests.unit.test_doc_reference_exit_codes import CHECKER, run_checker
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    FREQUENCY_HZ,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

REPO_ROOT = Path(__file__).resolve().parents[2]
PAGE = REPO_ROOT / "docs" / "guide" / "first_study.md"
TAG = "first-study"
OUTPUT_TAG = "first-study-output"

# `WF-15`'s degree-2 Study record (§9 item 251): the `.10e` print at
# `20261005T170537Z_WF-15-r2-W2a-deg2-study.log:170`; full repr
# 8.011657724511505e-08 at `:163`. Not imported: the WF-15 gate module holds
# no named constant for it, and adding one would force a re-run of its W2a.
WF15_DEG2_RECORD = 8.0116577245e-08

LINK_PASS_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")
OUTPUT_LINE_RE = re.compile(r"^mean SAR = (\S+) W/kg$")


def extract_tagged_block(text: str, tag: str, lang: str = "python") -> str:
    """The single fenced block opened by exactly "```<lang> <tag>"."""
    opener = f"```{lang} {tag}"
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line == opener]
    if len(starts) != 1:
        raise LookupError(f"{len(starts)} blocks tagged {opener!r}; need exactly one")
    body = []
    for line in lines[starts[0] + 1 :]:
        if line == "```":
            return "\n".join(body) + "\n"
        body.append(line)
    raise LookupError(f"block tagged {opener!r} is never closed")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _page_and_block():
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), TAG)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[WF-17] page sha256  = {_sha(raw)}")
        print(f"[WF-17] block sha256 = {_sha(block.encode())}")
    return raw.decode(), block


def _study_calls(source: str):
    calls = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call):
            fn = node.func
            name = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)
            if name == "Study":
                calls.append(node)
    return calls


class _StopBeforeSolve(Exception):
    pass


def exec_block(source: str, monkeypatch, *, stop_before_solve: bool):
    """Exec ``source`` as __main__ on every rank; return (recorded, stdout).

    A spy wraps ``Study.mean_sar`` (calls the original, records instance and
    dict); with ``stop_before_solve`` ``Study.solve`` records ``self`` and
    raises. The spy never alters what the block computes.
    """
    recorded: dict = {}
    original_mean_sar = Study.mean_sar

    def spy_mean_sar(self):
        recorded["study"] = self
        out = original_mean_sar(self)
        recorded["sar"] = out
        return out

    monkeypatch.setattr(Study, "mean_sar", spy_mean_sar)
    if stop_before_solve:

        def stop_solve(self, comm=MPI.COMM_WORLD):
            recorded["study"] = self
            raise _StopBeforeSolve

        monkeypatch.setattr(Study, "solve", stop_solve)

    buf = io.StringIO()
    code = compile(source, "docs/guide/first_study.md#first-study", "exec")
    with contextlib.redirect_stdout(buf):
        try:
            exec(code, {"__name__": "__main__"})
        except _StopBeforeSolve:
            if not stop_before_solve:
                raise
    captured = buf.getvalue()
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-17] block stdout (rank 0): {captured!r}")
    return recorded, captured


def configuration_mismatches(study) -> list[str]:
    """(K)'s comparator: which of WF-15's Study fields differ, exactly."""
    expected = {
        "geometry": GEOMETRY,
        "material": GelledSalinePhantomMaterial(SIGMA_HIGH, EPSILON_R_SPHERE, FREQUENCY_HZ),
        "excitation": QuasiStaticSphereField(e0=E0),
        "frequency_hz": FREQUENCY_HZ,
        "density_kg_m3": RHO_KG_M3,
        "degree": 2,
    }
    return [name for name, value in expected.items() if getattr(study, name) != value]


def _assert_extraction(block):
    assert block.strip(), "(X) the first-study block is empty"
    calls = _study_calls(block)
    assert len(calls) == 1, f"(X) {len(calls)} Study( calls in the block, need exactly 1"
    return calls[0]


# ---------------------------------------------------------------- W0 (smoke)


def test_wf17_extracts_exactly_one_first_study_block():
    _, block = _page_and_block()
    _assert_extraction(block)


def test_wf17_block_configures_the_wf15_study(monkeypatch):
    _, block = _page_and_block()
    call = _assert_extraction(block)
    keywords = [kw.arg for kw in call.keywords]
    print(f"[WF-17] Study keywords in the block: {keywords}")
    assert "degree" not in keywords, "(K) the tutorial must teach the default degree"
    recorded, _ = exec_block(block, monkeypatch, stop_before_solve=True)
    assert "study" in recorded, "(K) the block never reached Study.solve"
    mismatches = configuration_mismatches(recorded["study"])
    print(f"[WF-17] (K) mismatches: {mismatches}")
    assert mismatches == []


def test_wf17_planted_page_and_parameter_controls(monkeypatch):
    text, block = _page_and_block()
    # (c1) zero tags and two tags both raise.
    untagged = text.replace(f"```python {TAG}\n", "```python\n")
    with pytest.raises(LookupError):
        extract_tagged_block(untagged, TAG)
    doubled = text + f"\n```python {TAG}\n{block}```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    # (c2) sigma 0.57 -> 0.75 is caught before any solve, and only as material.
    assert block.count("0.57") == 1, "(c2) the sigma literal must occur exactly once"
    planted = block.replace("0.57", "0.75")
    recorded, _ = exec_block(planted, monkeypatch, stop_before_solve=True)
    mismatches = configuration_mismatches(recorded["study"])
    print(f"[WF-17] (c2) mismatches on sigma 0.75: {mismatches}")
    assert mismatches == ["material"]


def test_wf17_expected_output_states_the_record():
    text, _ = _page_and_block()
    out = extract_tagged_block(text, OUTPUT_TAG, lang="text")
    values = [OUTPUT_LINE_RE.match(line) for line in out.splitlines()]
    values = [m.group(1) for m in values if m]
    print(f"[WF-17] (E) expected-output values: {values}")
    assert len(values) == 1
    assert float(values[0]) == WF15_DEG2_RECORD


def _checker_stdout(*extra):
    proc = subprocess.run(
        [sys.executable, str(CHECKER), *extra], capture_output=True, text=True, cwd=REPO_ROOT
    )
    return proc.stdout


def test_wf17_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"guide={counts['guide']} stale={counts['stale']} exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(1)) >= 1


def test_wf17_docrefs_detects_a_planted_dead_link(tmp_path):
    if MPI.COMM_WORLD.rank != 0:
        return
    guides = tmp_path / "guides"
    guides.mkdir()
    (guides / "sibling.md").write_text("# sibling\n")
    (guides / "p.md").write_text(
        "# planted\n\n"
        "[live](sibling.md)\n\n"
        "[x](no_such_page.md)\n\n"
        "```text\n[y](also_missing.md)\n```\n"
    )
    _, with_links = run_checker("--docs-root", str(guides), "--check-links")
    _, without = run_checker("--docs-root", str(guides))
    print(f"[WF-17] (c3) with --check-links dead={with_links['dead']}; without dead={without['dead']} (printed)")
    assert with_links["dead"] == 1


# ------------------------------------------------------- complex: W1 and W2


@complex_only
def test_wf17_planted_degree1_block_misses_the_record(monkeypatch):
    _, block = _page_and_block()
    tree = ast.parse(block)
    calls = _study_calls(block)
    assert len(calls) == 1
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "Study":
            node.keywords.append(ast.keyword(arg="degree", value=ast.Constant(1)))
    planted = ast.unparse(ast.fix_missing_locations(tree))
    t0 = time.perf_counter()
    recorded, _ = exec_block(planted, monkeypatch, stop_before_solve=False)
    seconds = time.perf_counter() - t0
    sar = float(recorded["sar"]["mean_sar_w_per_kg"])
    assert recorded["study"].degree == 1
    rel_record = abs(sar - WF15_DEG2_RECORD) / WF15_DEG2_RECORD
    rel_mat4 = abs(sar - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[WF-17] (c4) degree-1 SAR = {sar!r} ({sar:.10e}); vs WF15 deg-2 record "
            f"{rel_record:.3e} (must be > {L0_BAND:g}); vs MAT4 deg-1 record "
            f"{rel_mat4:.3e} (must be <= {L0_BAND:g}); exec {seconds:.1f} s"
        )
    assert rel_record > L0_BAND
    assert rel_mat4 <= L0_BAND


@complex_only
def test_wf17_first_study_block_reproduces_wf15_record(monkeypatch):
    comm = MPI.COMM_WORLD
    _, block = _page_and_block()
    call = _assert_extraction(block)
    assert "degree" not in [kw.arg for kw in call.keywords]
    t0 = time.perf_counter()
    recorded, captured = exec_block(block, monkeypatch, stop_before_solve=False)
    seconds = time.perf_counter() - t0
    study = recorded["study"]
    # (K) on the live, unstubbed instance.
    assert configuration_mismatches(study) == []
    # (Deg)
    space = study.result.fields.e_complex.function_space
    element_degree = space.ufl_element().degree
    element_degree = element_degree() if callable(element_degree) else element_degree
    cells = study.result.mesh.topology.index_map(3).size_global
    dofs = space.dofmap.index_map.size_global * space.dofmap.index_map_bs
    if comm.rank == 0:
        print(
            f"[WF-17] (Deg) study.degree={study.degree} result.degree={study.result.degree} "
            f"element degree={element_degree}; cells={cells} N1curl dofs={dofs} "
            f"(printed; version-tagged 74020 / 482158)"
        )
    assert study.degree == 2
    assert study.result.degree == 2
    assert element_degree == 2
    # (R)
    sar = float(recorded["sar"]["mean_sar_w_per_kg"])
    rel_record = abs(sar - WF15_DEG2_RECORD) / WF15_DEG2_RECORD
    # (G)
    closed = _closed(SIGMA_HIGH, EPSILON_R_SPHERE)["sar_w_per_kg"]
    rel_closed = abs(sar - closed) / closed
    if comm.rank == 0:
        print(
            f"[WF-17] (R) SAR = {sar!r} ({sar:.10e}); record {WF15_DEG2_RECORD:.10e}; "
            f"rel {rel_record:.3e} (band L0_BAND {L0_BAND:g})"
        )
        print(
            f"[WF-17] (G) closed form {closed:.6e}; rel {rel_closed:.4%} (band MAT4_BAND {MAT4_BAND:g})"
        )
        print(f"[WF-17] block exec incl. mesh + degree-2 solve + SAR: {seconds:.1f} s")
    assert rel_record <= L0_BAND
    assert rel_closed <= MAT4_BAND
    # (Out)
    if comm.rank == 0:
        assert f"{sar:.10e}" in captured
