"""`WF-17` step 12: the copper-from-YAML tutorial page runs `WF-1`'s
f_small_copper.yaml through ``run_study_file`` and reproduces `ANS-6`'s
degree-1 copper class means at their band.

Brief: docs/planning/chunks/WF-17.md, "§9 item 619".  The solve node runs only
with ``FEM_EM_WF17_S12=1`` (and then asserts the complex build and ``-n
ANS6_EXACT_STRIP_RECORD_RANK_WIDTH``, never skips).  The other nodes build no
mesh.

Anchors: (X) one block per tag, one ``run_study_file(`` call, no ``Study(``,
the file argument read from the AST, the bash line's tokens, the YAML's
``conductor`` mapping = ``{"kind": "leontovich", "sigma_s_per_m":
COPPER_SIGMA}`` by ``float.hex``; (E) the bands block = the four imported bands
at ``.0e``; (R-run) one call, family / degree / summary conductor / Touchstone
name; (R) cells = ``ANS6_EXACT_STRIP_HOLE_CELL_RECORD`` and the class means
within ``LEG_D0_REPRODUCTION_BAND`` of the copper record; (g) TH-14's
reciprocity / passivity / spread gates on the run's S; (D) 0 dead links.
Controls (rule (e)): (c1) the PEC-control file (run by this module, never on
the page) ties the ``"pec"`` record and misses the copper record by more than
``NEGATIVE_CONTROL_FLOOR``, asserted, backed; (c-page) by construction.  Bands
only are asserted; every reading is printed, never asserted (item 574's
ruling).
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

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver
from fem_em_solver.study import COPPER_SIGMA

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_study_copper import NEGATIVE_CONTROL_FLOOR, _rel_means
from tests.package.test_wf17_first_study_tutorial import LINK_PASS_RE, _checker_stdout, extract_tagged_block
from tests.package.test_wf1_study_yaml import _mapping
from tests.package.test_wf1_yaml_copper import COPPER_FILE, KEY, PEC_CONTROL_FILE
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "copper_birdcage_from_yaml.md"
PAGE_ID = "docs/guide/copper_birdcage_from_yaml.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s12"
S12_ENV = "FEM_EM_WF17_S12"
MIN_GUIDE_PAGES = 61

COPPER_TAG = "copper-yaml"
CLI_TAG = "copper-yaml-cli"
BANDS_TAG = "copper-yaml-bands"
CLI_TOKENS = ["femem", "run", "f_small_copper.yaml", "--output-dir", "femem_out/f_small_copper"]

# Printed beside the readings (`20261020T202415Z_WF-1-s21-RGCu-copper.log`
# :2714 / :2716 / :2718); never asserted.
RGCU_R = "4.093e-11 / 8.673e-11 / 1.376e-10"
RGCU_C1 = "6.422e-04 / 2.812e-04 / 4.863e-04"
RGCU_G = "reciprocity 6.478e-16, max sigma 0.999851545099"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _calls_named(block: str, name: str) -> list:
    out = []
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Call):
            fn = node.func
            fname = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)
            if fname == name:
                out.append(node)
    return out


def _file_arg(block: str) -> str:
    """The string argument of the block's single ``run_study_file`` call (AST, never exec)."""
    calls = _calls_named(block, "run_study_file")
    if len(calls) != 1:
        raise LookupError(f"{len(calls)} run_study_file calls; need exactly one")
    return ast.literal_eval(calls[0].args[0])


def _bands_lines() -> list:
    return [
        f"S class means vs the copper record: {LEG_D0_REPRODUCTION_BAND:.0e}",
        f"Reciprocity, relative: {RECIPROCITY_BAND:.0e}",
        f"Largest singular value above 1: {PASSIVITY_SIGMA_TOLERANCE:.0e}",
        f"Class spread: {ADJACENT_SPREAD_BAND:.0e}",
    ]


def _bands_defect(block: str):
    """None when the bands block equals the four rendered lines, else the reason."""
    got = block.strip("\n").splitlines()
    want = _bands_lines()
    if len(got) != len(want):
        return f"{len(got)} lines != {len(want)}"
    for g, w in zip(got, want):
        if g != w:
            return f"line {g!r} != {w!r}"
    return None


def _knob() -> bool:
    raw = os.environ.get(S12_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S12_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(f"{S12_ENV} unset: the tutorial solve window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the copper four-port port Study needs the complex build"
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def _fmt(means):
    return ", ".join(f"{c} {means[c].real:+.12e}{means[c].imag:+.12e}j" for c in CLASSES)


def exec_copper_block(tag: str, monkeypatch):
    """Copy the block's file into the scratch, chdir there, exec the block.

    Returns ``(namespace, calls, stdout, scratch)``; the spy calls through.
    """
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), tag)
    name = _file_arg(block)
    scratch = SCRATCH_ROOT / tag
    if comm.rank == 0:
        print(f"[WF-17 s12] page sha256 {_sha(raw)}; block {tag} sha256 {_sha(block.encode())}",
              flush=True)
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        data = (COPPER_FILE.parent / name).read_bytes()
        (scratch / name).write_bytes(data)
        print(f"[WF-17 s12] copied {name} sha256 {_sha(data)} (tracked "
              f"{_sha(COPPER_FILE.read_bytes())})", flush=True)
    comm.Barrier()
    monkeypatch.chdir(scratch)
    calls = {"n": 0, "paths": []}
    original = fem_em_solver.run_study_file

    def spy(path, *args, **kwargs):
        calls["n"] += 1
        calls["paths"].append(str(path))
        return original(path, *args, **kwargs)

    monkeypatch.setattr(fem_em_solver, "run_study_file", spy)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{tag}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s12] block {tag} stdout (rank 0): {out!r}")
    return namespace, calls, out, scratch


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s12_blocks_and_file():
    text = _text()
    block = extract_tagged_block(text, COPPER_TAG)
    n_run = len(_calls_named(block, "run_study_file"))
    n_study = len(_calls_named(block, "Study"))
    name = _file_arg(block)
    _p0(f"[WF-17 s12] (X) {COPPER_TAG}: run_study_file calls {n_run}, Study calls {n_study}; "
        f"file {name!r} (want {COPPER_FILE.name!r})")
    assert block.count("run_study_file(") == 1
    assert n_run == 1 and n_study == 0 and "Study(" not in block
    assert name == COPPER_FILE.name
    assert COPPER_FILE.is_file() and PEC_CONTROL_FILE.is_file()
    assert f"(../../tests/package/wf1/{COPPER_FILE.name})" in text
    cli = extract_tagged_block(text, CLI_TAG, lang="bash")
    cli_lines = [ln for ln in cli.splitlines() if ln.strip()]
    _p0(f"[WF-17 s12] (X) {CLI_TAG}: {cli_lines!r}")
    assert len(cli_lines) == 1
    assert cli_lines[0].split() == CLI_TOKENS
    extract_tagged_block(text, BANDS_TAG, lang="text")
    # the YAML's conductor mapping (parse only, no mesh)
    cond = _mapping(COPPER_FILE)["conductor"]
    _p0(f"[WF-17 s12] (X) conductor {cond!r}; sigma hex "
        f"{float(cond['sigma_s_per_m']).hex()} (COPPER_SIGMA {COPPER_SIGMA.hex()})")
    assert set(cond) == {"kind", "sigma_s_per_m"}
    assert cond["kind"] == "leontovich"
    assert float(cond["sigma_s_per_m"]).hex() == COPPER_SIGMA.hex()
    assert "conductor" not in _mapping(PEC_CONTROL_FILE)
    # (c-page), by construction: a removed / doubled tag raises
    for tag, lang in ((COPPER_TAG, "python"), (CLI_TAG, "bash"), (BANDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    _p0("[WF-17 s12] (c-page) removed / doubled tags raise LookupError for all three blocks")
    for token in ("Study(", ".json", ".npz", ".s4p", ".xdmf"):
        assert token not in text.replace("run_study_file(", ""), token


def test_wf17_s12_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s12] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _bands_defect(block) is None, _bands_defect(block)
    # (c-page), by construction: a line at a different exponent fails (E) and is named
    altered = block.replace(f"{RECIPROCITY_BAND:.0e}", f"{RECIPROCITY_BAND * 10:.0e}", 1)
    reason = _bands_defect(altered)
    _p0(f"[WF-17 s12] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "Reciprocity, relative" in reason


def test_wf17_s12_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s12] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s12_copper_block_reproduces_ans6(monkeypatch):
    comm = _solve_preamble()
    assert th14.FREQUENCIES_HZ[KEY] == FREQUENCY_64_HZ
    t0 = time.perf_counter()
    ns, calls, _, scratch = exec_copper_block(COPPER_TAG, monkeypatch)
    t_block = time.perf_counter() - t0
    run = ns["run"]
    s = run.summary

    # (R) the copper record through the page
    means = _class_means(run.result.s_matrix)
    rel = _rel_means(means, ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    # (g) TH-14's gates on the page's S
    rec = th14._sweep_record(run.result.sweep)
    sigma_max = float(np.max(rec["sigma"]))

    # (c1) the PEC control, run here (never on the page), output under the scratch
    n_calls, paths = calls["n"], list(calls["paths"])
    t1 = time.perf_counter()
    run_pec = fem_em_solver.run_study_file(PEC_CONTROL_FILE, output_dir=scratch / "pec", comm=comm)
    t_pec = time.perf_counter() - t1
    pec = _class_means(run_pec.result.s_matrix)
    rel_pec_tie = _rel_means(pec, ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["pec"])
    rel_c1 = _rel_means(pec, ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])

    _p0(f"[WF-17 s12] (R-run) block spy calls {n_calls} {paths}; family {s['family']!r}; "
        f"degree {s['degree']!r}; conductor {s.get('conductor', '<absent>')!r}; touchstone "
        f"{run.paths['touchstone'].name!r}")
    _p0(f"[WF-17 s12] (R) cells {run.result.cells} (record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); "
        f"cu class means {_fmt(means)}")
    _p0("[WF-17 s12] (R) page cu vs record 'cu', relative: " + " / ".join(
        f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}; digits printed, never asserted; "
        f"RGCu:2714 {RGCU_R})")
    _p0(f"[WF-17 s12] (g) reciprocity {rec['reciprocity']:.3e} (ASSERTED <= "
        f"{RECIPROCITY_BAND:.0e}); sigma {[f'{v:.12f}' for v in rec['sigma']]}, max "
        f"{sigma_max:.12f} (ASSERTED <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e}); spreads "
        f"{rec['spreads']} (ASSERTED <= {ADJACENT_SPREAD_BAND:.0e}); RGCu:2718 {RGCU_G}")
    _p0(f"[WF-17 s12] (c1) PEC control cells {run_pec.result.cells}; class means {_fmt(pec)}")
    _p0("[WF-17 s12] (c1) PEC control vs record 'pec' (tie, printed): " + " / ".join(
        f"{c} {rel_pec_tie[c]:.3e}" for c in CLASSES)
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})")
    _p0("[WF-17 s12] (c1) PEC control vs record 'cu' (miss): " + " / ".join(
        f"{c} {rel_c1[c]:.3e}" for c in CLASSES)
        + f" (ASSERTED > {NEGATIVE_CONTROL_FLOOR:.0e}; RGCu:2716 {RGCU_C1}); min/floor "
        f"{min(rel_c1.values()) / NEGATIVE_CONTROL_FLOOR:.1f}x")
    _p0(f"[WF-17 s12] W-cu block (cu run_study_file) {t_block:.1f} s; PEC control "
        f"run_study_file {t_pec:.1f} s; node {time.perf_counter() - t0:.1f} s")

    # (R-run)
    assert n_calls == 1
    assert paths == [COPPER_FILE.name]
    assert s["family"] == "four_port_birdcage"
    assert s["degree"] == 1
    assert s["conductor"] == {"kind": "leontovich", "sigma_s_per_m": COPPER_SIGMA}
    assert run.paths["touchstone"].name == COPPER_FILE.stem + ".s4p"
    assert "conductor" not in run_pec.summary
    # (R)
    assert run.result.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])
    # (g)
    assert rec["reciprocity"] <= RECIPROCITY_BAND, rec["reciprocity"]
    assert sigma_max <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, rec["sigma"]
    for name, value in rec["spreads"].items():
        assert value <= ADJACENT_SPREAD_BAND, (name, value)
    # (c1)
    for c in CLASSES:
        assert rel_pec_tie[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel_pec_tie[c])
        assert rel_c1[c] > NEGATIVE_CONTROL_FLOOR, (c, rel_c1[c])
