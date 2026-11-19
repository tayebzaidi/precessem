"""`WF-17` step 21 (§9 item 654, U93): the Cole–Cole sphere page is a test.

The committed page `docs/guide/dispersive_tissue_sphere.md` carries one fenced
block tagged ```` ```python dispersive-sphere ````. This module reads the
committed page, extracts that block verbatim and ``exec``s it as ``__main__``
(degree 2 by default, ``-n 2``, complex). Brief:
docs/planning/chunks/WF-17.md "§9 item 654". Rule (i): this module and the
page are not edited after their last green window.

Asserted:
(X)    exactly one tagged block with exactly one ``Study(`` call and no
       ``degree=`` keyword; exec'd at -n 2 it reproduces U92's W1 SAR record at
       every digit that record carries (see RECORD below for why that is the
       bound, not ``L0_BAND``);
(Med)  the page's medium equals `MAT-3` step 1b's ``COLE_COLE``: ε∞, the poles
       and σ_i exact (float equality), and the Study is U92's otherwise;
(E)    the page's expected-output block states U92's record and the resolved
       (σ_eff, εᵣ) as the block prints them;
(G)    the page's SAR is within ``MAT4_BAND`` of ``_closed(σ_eff, εᵣ)`` at
       ε̂(64 MHz);
(D)    ``docrefs --docs-root docs/guide --check-links`` reads dead == 0.

Negative controls (asserted): (c-page) a tag removed or doubled raises; a
planted σ_i 0.2 -> 0.3 fails (Med) and only as the medium; the page's SAR
misses the Debye-reduced closed form by > ``MAT4_BAND`` (U92's (c-D)); a
planted dead link reads dead == 1.
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

from fem_em_solver import QuasiStaticSphereField, Study
from tests.complex_mode import complex_only
from tests.materials.test_phantom_material_model import MAT4_BAND
from tests.package.test_wf15_study_lossy_sphere import GEOMETRY, _closed
from tests.unit.test_doc_reference_exit_codes import CHECKER, run_checker
from tests.validation.test_lossy_sphere_sar import E0, FREQUENCY_HZ, RHO_KG_M3
from tests.validation.test_mat3_cole_cole_plane_wave import (
    COLE_COLE,
    DEBYE_REDUCED,
    EPS_INF,
    POLES,
    SIGMA_IONIC,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

REPO_ROOT = Path(__file__).resolve().parents[2]
PAGE = REPO_ROOT / "docs" / "guide" / "dispersive_tissue_sphere.md"
TAG = "dispersive-sphere"
OUTPUT_TAG = "dispersive-sphere-output"

# U92's (`WF-15` step 8, §9 item 653) degree-2 Cole–Cole Study SAR, by value:
# `20261027T003215Z_WF-15-s8-W1.log:135` prints "FEM = 3.986757e-08" (.6e).
# That print is the only record: the gate holds no named constant and printed
# no fuller repr. The brief asks for L0_BAND (1e-9 relative), but a 7-digit
# record carries a rounding half-width of 0.5e-14 / 3.986757e-08 = 1.254e-07
# relative, 125x wider than L0_BAND, so L0_BAND against this number would test
# its rounding, not the page. The bound used is the record's own precision:
# the page's SAR formatted at .6e equals RECORD_STR, digit for digit (the
# "digits equal" form of steps 8 / 9 / 12). The L0_BAND distance is printed.
RECORD_STR = "3.986757e-08"
RECORD = float(RECORD_STR)

LINK_PASS_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")
SAR_LINE_RE = re.compile(r"^mean SAR = (\S+) W/kg$")


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


def _page_and_block():
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), TAG)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[WF-17 s21] page sha256  = {hashlib.sha256(raw).hexdigest()}")
        print(f"[WF-17 s21] block sha256 = {hashlib.sha256(block.encode()).hexdigest()}")
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


def _assert_extraction(block):
    assert block.strip(), "(X) the dispersive-sphere block is empty"
    calls = _study_calls(block)
    assert len(calls) == 1, f"(X) {len(calls)} Study( calls in the block, need exactly 1"
    assert "degree" not in [
        kw.arg for kw in calls[0].keywords
    ], "(X) the page teaches the default degree"
    return calls[0]


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

        def stop_solve(self, comm=MPI.COMM_WORLD, **kwargs):
            recorded["study"] = self
            raise _StopBeforeSolve

        monkeypatch.setattr(Study, "solve", stop_solve)

    buf = io.StringIO()
    code = compile(source, "docs/guide/dispersive_tissue_sphere.md#dispersive-sphere", "exec")
    with contextlib.redirect_stdout(buf):
        try:
            exec(code, {"__name__": "__main__"})
        except _StopBeforeSolve:
            if not stop_before_solve:
                raise
    captured = buf.getvalue()
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-17 s21] block stdout (rank 0): {captured!r}")
    return recorded, captured


def medium_mismatches(study) -> list[str]:
    """(Med)'s comparator: which of U92's Study fields differ, exactly."""
    med = study.dispersive_material
    out = []
    if med is None:
        return ["medium"]
    if med.epsilon_inf != EPS_INF or med.poles != POLES or med.sigma_ionic != SIGMA_IONIC:
        out.append("medium")
    if med != COLE_COLE:
        out.append("medium-dataclass")
    expected = {
        "geometry": GEOMETRY,
        "excitation": QuasiStaticSphereField(e0=E0),
        "frequency_hz": FREQUENCY_HZ,
        "density_kg_m3": RHO_KG_M3,
        "degree": 2,
    }
    out += [name for name, value in expected.items() if getattr(study, name) != value]
    return out


def _resolved_line(sigma: float, eps_r: float) -> str:
    return f"resolved at 64 MHz: sigma = {sigma:.6g} S/m, epsilon_r = {eps_r:.6g}"


# ---------------------------------------------------------------- W0 (no solve)


def test_wf17_s21_block_medium_and_controls(monkeypatch):
    """(X)-extract, (Med) with the solve stopped, (c-page) tags and σ_i."""
    text, block = _page_and_block()
    _assert_extraction(block)
    recorded, captured = exec_block(block, monkeypatch, stop_before_solve=True)
    assert "study" in recorded, "(Med) the block never reached Study.solve"
    mism = medium_mismatches(recorded["study"])
    sigma_eff, eps_r = COLE_COLE.sigma_epsilon_at(FREQUENCY_HZ)
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[WF-17 s21] (Med) mismatches {mism}; medium {recorded['study'].dispersive_material!r}; "
            f"resolved sigma {recorded['study'].material.sigma!r} eps_r "
            f"{recorded['study'].material.epsilon_r!r}"
        )
    assert mism == []
    assert recorded["study"].material.sigma == sigma_eff
    assert recorded["study"].material.epsilon_r == eps_r
    assert captured.splitlines()[0] == _resolved_line(sigma_eff, eps_r)

    # (c-page) zero tags and two tags both raise.
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {TAG}\n", "```python\n"), TAG)
    with pytest.raises(LookupError):
        extract_tagged_block(text + f"\n```python {TAG}\n{block}```\n", TAG)
    with pytest.raises(LookupError):
        extract_tagged_block(
            text.replace(f"```text {OUTPUT_TAG}\n", "```text\n"), OUTPUT_TAG, "text"
        )
    # (c-page) sigma_ionic 0.2 -> 0.3 is caught as the medium only.
    assert block.count("sigma_ionic=0.2,") == 1
    planted = block.replace("sigma_ionic=0.2,", "sigma_ionic=0.3,")
    rec2, _ = exec_block(planted, monkeypatch, stop_before_solve=True)
    m2 = medium_mismatches(rec2["study"])
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-17 s21] (c-page) sigma_ionic 0.3 mismatches: {m2}")
    assert m2 == ["medium", "medium-dataclass"]


def test_wf17_s21_expected_output_states_the_record():
    text, _ = _page_and_block()
    out = extract_tagged_block(text, OUTPUT_TAG, lang="text").splitlines()
    sigma_eff, eps_r = COLE_COLE.sigma_epsilon_at(FREQUENCY_HZ)
    values = [m.group(1) for m in (SAR_LINE_RE.match(line) for line in out) if m]
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-17 s21] (E) output block {out}")
    assert out == [_resolved_line(sigma_eff, eps_r), f"mean SAR = {RECORD_STR} W/kg"]
    assert values == [RECORD_STR]


def _checker_stdout(*extra):
    proc = subprocess.run(
        [sys.executable, str(CHECKER), *extra], capture_output=True, text=True, cwd=REPO_ROOT
    )
    return proc.stdout


def test_wf17_s21_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s21] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(1)) >= 1


def test_wf17_s21_docrefs_detects_a_planted_dead_link(tmp_path):
    if MPI.COMM_WORLD.rank != 0:
        return
    guides = tmp_path / "guides"
    guides.mkdir()
    (guides / "sibling.md").write_text("# sibling\n")
    (guides / "p.md").write_text("# planted\n\n[live](sibling.md)\n\n[x](no_such_page.md)\n")
    _, with_links = run_checker("--docs-root", str(guides), "--check-links")
    print(f"[WF-17 s21] (c-dead) planted dead={with_links['dead']}")
    assert with_links["dead"] == 1


# ------------------------------------------------------- W-x (complex, -n 2)


@complex_only
def test_wf17_s21_dispersive_block_reproduces_u92_record(monkeypatch):
    comm = MPI.COMM_WORLD
    _, block = _page_and_block()
    _assert_extraction(block)
    t0 = time.perf_counter()
    recorded, captured = exec_block(block, monkeypatch, stop_before_solve=False)
    seconds = time.perf_counter() - t0
    study = recorded["study"]
    # (Med) on the live instance; degree 2.
    assert medium_mismatches(study) == []
    assert study.result.degree == 2
    cells = study.result.mesh.topology.index_map(3).size_global
    sar = float(recorded["sar"]["mean_sar_w_per_kg"])
    # (X)
    rel_record = abs(sar - RECORD) / RECORD
    # (G)
    sigma_eff, eps_r = COLE_COLE.sigma_epsilon_at(FREQUENCY_HZ)
    closed = _closed(sigma_eff, eps_r)["sar_w_per_kg"]
    rel_closed = abs(sar - closed) / closed
    # (c-page) Debye-reduced closed form.
    sd, ed = DEBYE_REDUCED.sigma_epsilon_at(FREQUENCY_HZ)
    closed_d = _closed(sd, ed)["sar_w_per_kg"]
    miss_d = abs(sar - closed_d) / closed_d
    if comm.rank == 0:
        print(
            f"[WF-17 s21] degree-{study.result.degree}, {cells} cells; sigma_eff {sigma_eff:.6g}, "
            f"eps_r {eps_r:.6g}\n"
            f"[WF-17 s21] (X) SAR = {sar!r} ({sar:.6e}); record {RECORD_STR} "
            f"(20261027T003215Z_WF-15-s8-W1.log:135); rel {rel_record:.3e} (printed vs "
            f"L0_BAND {L0_BAND:g}; asserted: .6e digits equal, half-width 1.254e-07)\n"
            f"[WF-17 s21] (G) closed {closed:.6e}; rel {rel_closed:.4%} (band MAT4_BAND {MAT4_BAND:g})\n"
            f"[WF-17 s21] (c-page) Debye-reduced closed {closed_d:.6e}; miss {miss_d:.4%} "
            f"(> {MAT4_BAND:g})\n"
            f"[WF-17 s21] block exec incl. mesh + degree-2 solve + SAR: {seconds:.1f} s"
        )
    assert f"{sar:.6e}" == RECORD_STR, f"(X) {sar:.6e} != record {RECORD_STR}"
    assert rel_closed <= MAT4_BAND, f"(G) {rel_closed:.2%}"
    assert miss_d > MAT4_BAND, f"(c-page) blind: within {miss_d:.2%} of Debye-reduced"
    # (Out) the block prints what it returns, as the output block states.
    assert captured.splitlines() == [
        _resolved_line(sigma_eff, eps_r),
        f"mean SAR = {sar:.6e} W/kg",
    ]
