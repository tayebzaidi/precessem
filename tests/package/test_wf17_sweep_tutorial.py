"""`WF-17` step 9: the sweep-from-YAML tutorial page runs `WF-1`'s
f_small_sweep.yaml through ``run_study_file`` and reproduces U71's three
per-frequency `ANS-4` P1 ties and the one-Touchstone round trip at their bands.

Brief: docs/planning/chunks/WF-17.md, "§9 item 591" and the review ruling
2026-10-22 03:00.  The solve node runs only with ``FEM_EM_WF17_S9=1`` (and then
asserts the complex build and ``-n ANS4_EXACT_STRIP_RECORD_RANK_WIDTH``, never
skips).  The other nodes build no mesh.

Anchors: (X) one block per tag, one ``run_study_file(`` call, no ``Study(``,
the file argument read from the AST, the YAML's ``sweep_hz`` = ``FREQUENCIES``
by ``float.hex``; (E) the bands block = the three imported bands at ``.0e``;
(R-run) one call, family ``four_port_birdcage``, no scalar solve; (R′) the
three ties at ``LEG_D0_REPRODUCTION_BAND`` on one mesh of
``ANS4_EXACT_STRIP_CELL_RECORD`` cells; (T′) the round trip at ``ANCHOR_TOL`` /
``FREQ_ATOL_HZ``; (D) 0 dead links.  Controls (rule (e)): (c-swap) asserted,
backed; (c-page) by construction.  Bands only are asserted; every reading is
printed, never asserted (item 574's ruling).
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
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study_yaml import _load

from tests.complex_mode import IS_COMPLEX
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ
from tests.package.test_wf16_study_four_port import _rel
from tests.package.test_wf16_study_sweep import C4_HAND_TABLE, KEYS, _records
from tests.package.test_wf17_first_study_tutorial import LINK_PASS_RE, _checker_stdout, extract_tagged_block
from tests.package.test_wf1_yaml_sweep import SWEEP, _unpairs
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    FREQUENCIES,
)
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "sweep_from_yaml.md"
PAGE_ID = "docs/guide/sweep_from_yaml.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s9"
S9_ENV = "FEM_EM_WF17_S9"
MIN_GUIDE_PAGES = 58

SWEEP_TAG = "four-port-sweep"
CLI_TAG = "four-port-sweep-cli"
BANDS_TAG = "four-port-sweep-bands"

# Printed beside the readings (U71 W1a `20261020T201454Z_WF-1-s21-W1a-sweep-file.log:1934/:1936`,
# EX-139 `20261021T020252Z_EX-139.log:1909/:1910`); never asserted.
W1A_R = "6.601e-11 / 9.557e-11 / 9.449e-11"
EX139_R = "6.601e-11 / 9.557e-11 / 9.449e-11"
W1A_T = "5.541e-14 / 6.873e-14 / 4.969e-14"
EX139_T = "5.456e-14 / 6.855e-14 / 4.969e-14"


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
        f"S vs the P1 records, per frequency: {LEG_D0_REPRODUCTION_BAND:.0e}",
        f"Touchstone round trip, S: {ANCHOR_TOL:.0e}",
        f"Touchstone round trip, frequency (Hz): {FREQ_ATOL_HZ:.0e}",
    ]


def _bands_defect(block: str):
    """None when the bands block equals the three rendered lines, else the reason."""
    got = block.strip("\n").splitlines()
    want = _bands_lines()
    if len(got) != len(want):
        return f"{len(got)} lines != {len(want)}"
    for g, w in zip(got, want):
        if g != w:
            return f"line {g!r} != {w!r}"
    return None


def _knob() -> bool:
    raw = os.environ.get(S9_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S9_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(f"{S9_ENV} unset: the tutorial solve window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port port Study needs the complex build"
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def exec_sweep_block(tag: str, monkeypatch):
    """Copy the block's file into the scratch, chdir there, exec the block.

    Returns ``(namespace, calls, stdout)``; the spy calls through.
    """
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), tag)
    name = _file_arg(block)
    scratch = SCRATCH_ROOT / tag
    if comm.rank == 0:
        print(f"[WF-17 s9] page sha256 {_sha(raw)}; block {tag} sha256 {_sha(block.encode())}",
              flush=True)
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        data = (SWEEP.parent / name).read_bytes()
        (scratch / name).write_bytes(data)
        print(f"[WF-17 s9] copied {name} sha256 {_sha(data)}", flush=True)
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
    _p0(f"[WF-17 s9] block {tag} stdout (rank 0): {out!r}")
    return namespace, calls, out


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s9_blocks_and_file():
    text = _text()
    block = extract_tagged_block(text, SWEEP_TAG)
    n_run = len(_calls_named(block, "run_study_file"))
    n_study = len(_calls_named(block, "Study"))
    name = _file_arg(block)
    _p0(f"[WF-17 s9] (X) {SWEEP_TAG}: run_study_file calls {n_run}, Study calls {n_study}; "
        f"file {name!r} (want {SWEEP.name!r})")
    assert block.count("run_study_file(") == 1
    assert n_run == 1 and n_study == 0 and "Study(" not in block
    assert name == SWEEP.name
    assert SWEEP.is_file()
    assert f"(../../tests/package/wf1/{SWEEP.name})" in text
    cli = extract_tagged_block(text, CLI_TAG, lang="bash")
    cli_lines = [ln for ln in cli.splitlines() if ln.strip()]
    _p0(f"[WF-17 s9] (X) {CLI_TAG}: {cli_lines!r}")
    assert len(cli_lines) == 1
    assert cli_lines[0].startswith("femem run ")
    assert SWEEP.name in cli_lines[0].split()
    extract_tagged_block(text, BANDS_TAG, lang="text")
    # the YAML's sweep_hz = FREQUENCIES by float.hex (parse only, no mesh)
    _study, kind, outputs, _cooling = _load(SWEEP)
    got = [float(f).hex() for f in outputs["sweep_hz"]]
    want = [float(f).hex() for f in FREQUENCIES.values()]
    _p0(f"[WF-17 s9] (X) kind {kind!r}; sweep_hz {outputs['sweep_hz']}; hex {got} (want {want})")
    assert kind == "four_port_birdcage"
    assert got == want
    # (c-page), by construction: a removed / doubled tag raises
    for tag, lang in ((SWEEP_TAG, "python"), (CLI_TAG, "bash"), (BANDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    _p0("[WF-17 s9] (c-page) removed / doubled tags raise LookupError for all three blocks")
    for token in ("Study(", ".json", ".npz", ".s4p", ".xdmf", "femem check"):
        assert token not in text.replace("run_study_file(", ""), token


def test_wf17_s9_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s9] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _bands_defect(block) is None, _bands_defect(block)
    # (c-page), by construction: a line at a different exponent fails (E) and is named
    altered = block.replace(f"{ANCHOR_TOL:.0e}", f"{ANCHOR_TOL * 10:.0e}", 1)
    reason = _bands_defect(altered)
    _p0(f"[WF-17 s9] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "Touchstone round trip, S" in reason


def test_wf17_s9_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s9] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s9_sweep_block_reproduces_u71(monkeypatch):
    _solve_preamble()
    t0 = time.perf_counter()
    ns, calls, _ = exec_sweep_block(SWEEP_TAG, monkeypatch)
    t_block = time.perf_counter() - t0
    run = ns["run"]
    s = run.summary
    sw = run.result
    _p0(f"[WF-17 s9] (R-run) spy calls {calls['n']} {calls['paths']}; family "
        f"{s['family']!r}; study.result is None {run.study.result is None}")

    # (R′) per frequency, one mesh (the gate's computation)
    rec = _records(1)
    s_file = np.array([_unpairs(m) for m in s["s_matrices"]])
    rels = {k: _rel(s_file[i], rec[k]) for i, k in enumerate(KEYS)}
    one_mesh = all(r.mesh is sw.results[0].mesh for r in sw.results)
    hex_ok = [float(f).hex() for f in s["sweep_hz"]] == [float(f).hex() for f in FREQUENCIES.values()]

    # (c-swap) asserted: the 10 <-> 128 MHz exchange misses
    i10, i128 = KEYS.index("10"), KEYS.index("128")
    sw_10 = _rel(s_file[i10], rec["128"])
    sw_128 = _rel(s_file[i128], rec["10"])

    # (T′) the one Touchstone
    ts = run.paths["touchstone"]
    freqs, mats, _ports = load_touchstone(ts)
    loaded = np.asarray(mats, dtype=np.complex128)
    n_pts = len(freqs)
    d_ts = [float(np.max(np.abs(loaded[i] - sw.s_matrices[i]))) for i in range(min(n_pts, len(KEYS)))]
    d_f = [abs(float(freqs[i]) - f) for i, f in enumerate(s["sweep_hz"][:n_pts])]

    _p0(f"[WF-17 s9] (R') ||S_f - P1_f||/||P1_f|| "
        + " / ".join(f"{k} MHz {v:.3e}" for k, v in rels.items())
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}; digits printed, never asserted; "
        f"W1a {W1A_R}; EX-139 {EX139_R}); cells {s['cells']} (record "
        f"{ANS4_EXACT_STRIP_CELL_RECORD}); one mesh {one_mesh}; sweep_hz hex equal {hex_ok}")
    _p0(f"[WF-17 s9] (c-swap) S_10 vs P1_128 {sw_10:.4f} (hand {C4_HAND_TABLE[1][('10', '128')]}); "
        f"S_128 vs P1_10 {sw_128:.4f} (hand {C4_HAND_TABLE[1][('128', '10')]}) "
        f"(ASSERTED > {LEG_D0_REPRODUCTION_BAND:.0e}; W1a 1.2367 / 1.2763)")
    _p0(f"[WF-17 s9] (T') {ts.name}: {n_pts} points; |df| {d_f} (ASSERTED <= {FREQ_ATOL_HZ:.0e}); "
        f"max|dS| " + " / ".join(f"{v:.3e}" for v in d_ts)
        + f" (ASSERTED <= {ANCHOR_TOL:.0e}; W1a {W1A_T}; EX-139 {EX139_T})")
    _p0(f"[WF-17 s9] W35 block {t_block:.1f} s (rank-0 summary elapsed {s['elapsed_s']:.1f} s; "
        f"mesh {s['mesh_time_s']:.2f} s; sweeps {s['sweep_times_s']}); node "
        f"{time.perf_counter() - t0:.1f} s")

    assert calls["n"] == 1
    assert calls["paths"] == [SWEEP.name]
    assert s["family"] == "four_port_birdcage"
    assert run.study.result is None
    assert s["cells"] == ANS4_EXACT_STRIP_CELL_RECORD
    assert one_mesh
    assert len(sw.results) == len(KEYS)
    assert hex_ok
    for k, v in rels.items():
        assert v <= LEG_D0_REPRODUCTION_BAND, k
    assert sw_10 > LEG_D0_REPRODUCTION_BAND
    assert sw_128 > LEG_D0_REPRODUCTION_BAND
    assert ts.name == SWEEP.stem + ".s4p" and s["touchstone"] == ts.name
    assert n_pts == len(KEYS) and loaded.shape == (len(KEYS), 4, 4)
    assert all(d <= FREQ_ATOL_HZ for d in d_f)
    assert all(d <= ANCHOR_TOL for d in d_ts)
