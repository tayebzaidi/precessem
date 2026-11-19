"""`WF-17` step 8: the four-port EM outputs tutorial page runs `WF-1`'s
f_small_sar_b1.yaml / f_small_em_probes.yaml through ``run_study_file`` and
reproduces item 344's 10 g / ``|B₁⁺|`` C4 readings and U70's identity.

Brief: docs/planning/chunks/WF-17.md, "§9 item 574" and the review ruling
2026-10-20 10:30.  The two solve nodes run only with ``FEM_EM_WF17_S8=1`` (and
then assert the complex build and ``-n ANS4_EXACT_STRIP_RECORD_RANK_WIDTH``,
never skip).  The other nodes build no mesh.

Anchors: (X) one block per tag, one ``run_study_file(`` call, no ``Study(``,
the file argument read from the AST; (M) the malformed files refuse naming
their field; (E) the records block = item 344's printed digits; (R-run) one
call per block, family ``four_port_birdcage``; (SB) the S tie at
``LEG_D0_REPRODUCTION_BAND``, the 10 g cyclic pairs and the ``|B₁⁺|``
covariance at ``C4_COVARIANCE_BAND`` and at their ``.4f`` digits; (PR) the
probe summary = ``point_sar`` on the run's fields on the same five-point
batch, the ``|E|`` identity at ``EXACT_IDENTITY_RTOL``; (D) 0 dead links.
Controls (rule (e)): (c1) and (c-disp) asserted, backed; (c-page) by
construction.  Every band and helper is imported, none restated; the seven
record digits are the printed readings (no named constant holds them).
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
from fem_em_solver import cli, load_study
from fem_em_solver.post.sar import point_sar
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_study_yaml import WF1
from tests.package.test_wf1_yaml_em_probes import EM_PROBES, POINTS, _field
from tests.package.test_wf1_yaml_sar_b1 import MALFORMED_SAR, SAR, SAR_B1_FILE, _diag
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import QUADRATURE_STEP_DEG, _port_index
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_post8_step3a_averaged_row1 import WIDTH_RTOL

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "four_port_outputs_from_yaml.md"
PAGE_ID = "docs/guide/four_port_outputs_from_yaml.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s8"
S8_ENV = "FEM_EM_WF17_S8"
MIN_GUIDE_PAGES = 57

SAR_TAG = "four-port-sar-b1"
PROBES_TAG = "four-port-probes"
CHECK_TAG = "four-port-check"
RECORDS_TAG = "four-port-records"
BLOCK_FILES = {SAR_TAG: SAR_B1_FILE.name, PROBES_TAG: EM_PROBES.name}

# Item 344's printed readings (W1 `20261009T005219Z_WF-1-s5-W1-sar-b1.log.gz`,
# re-read at equal digits by RG-SB `20261020T124017Z_WF-1-s20-RGSB-sar-b1.log`).
PAIRS_LABEL = "10 g C4 cyclic pairs (%)"
B1_LABEL = "|B1+| C4 covariance (%)"
PAIRS_RECORD = ["0.8702", "3.5284", "2.0337", "4.5024"]  # 0->1 / 1->2 / 2->3 / 3->0
B1_RECORD = ["1.8695", "1.9494", "1.6405"]  # P2@+D / P4@-D / P3@2D
DISPLACED_0 = (0.021, 0.0, 0.0)  # point 0 moved 1 mm radially (c-disp), as the gate


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


def _records(block: str) -> dict:
    out = {}
    for line in block.splitlines():
        if "): " in line:
            label, values = line.split("): ", 1)
            out[label + ")"] = values.split()
    return out


def _records_defect(block: str):
    """None when the records block states item 344's digits, else the reason."""
    got = _records(block)
    want = {PAIRS_LABEL: PAIRS_RECORD, B1_LABEL: B1_RECORD}
    if set(got) != set(want):
        return f"labels {sorted(got)} != {sorted(want)}"
    for label, values in want.items():
        if got[label] != values:
            return f"{label}: {got[label]} != {values}"
    return None


def _knob() -> bool:
    raw = os.environ.get(S8_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S8_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(
            f"{S8_ENV} unset: the tutorial solve windows are environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port port Study needs the complex build"
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def exec_four_port_block(tag: str, monkeypatch):
    """Copy the block's file into the scratch, chdir there, exec the block.

    Returns ``(namespace, calls, stdout)``; the spy calls through.
    """
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), tag)
    name = _file_arg(block)
    scratch = SCRATCH_ROOT / tag
    if comm.rank == 0:
        print(
            f"[WF-17 s8] page sha256 {_sha(raw)}; block {tag} sha256 {_sha(block.encode())}",
            flush=True,
        )
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        data = (WF1 / name).read_bytes()
        (scratch / name).write_bytes(data)
        print(f"[WF-17 s8] copied {name} sha256 {_sha(data)}", flush=True)
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
    _p0(f"[WF-17 s8] block {tag} stdout (rank 0): {out!r}")
    return namespace, calls, out


def _check_r_run(namespace, calls, tag):
    run = namespace["run"]
    _p0(
        f"[WF-17 s8] (R-run) {tag}: spy calls {calls['n']} {calls['paths']}; family "
        f"{run.summary['family']!r}"
    )
    assert calls["n"] == 1
    assert calls["paths"] == [BLOCK_FILES[tag]]
    assert run.summary["family"] == "four_port_birdcage"
    return run


def _s_tie(run):
    s = np.array([[complex(re, im) for re, im in row] for row in run.summary["s_matrix"]])
    rel = _rel(s, _record(1))
    _p0(
        f"[WF-17 s8] (R) ||S_yaml - P1|| / ||P1|| = {rel:.3e} (ASSERTED <= "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}; digits printed, never asserted)"
    )
    return rel


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s8_blocks_and_files():
    text = _text()
    for tag, want in BLOCK_FILES.items():
        block = extract_tagged_block(text, tag)
        n_run = len(_calls_named(block, "run_study_file"))
        n_study = len(_calls_named(block, "Study"))
        name = _file_arg(block)
        _p0(
            f"[WF-17 s8] (X) {tag}: run_study_file calls {n_run}, Study calls {n_study}; "
            f"file {name!r} (want {want!r})"
        )
        assert block.count("run_study_file(") == 1
        assert n_run == 1 and n_study == 0 and "Study(" not in block
        assert name == want
        assert (WF1 / name).is_file(), name
        assert f"(../../tests/package/wf1/{name})" in text, name
    extract_tagged_block(text, CHECK_TAG, lang="bash")
    extract_tagged_block(text, RECORDS_TAG, lang="text")
    # (c-page), by construction: a removed / doubled tag raises
    for tag in BLOCK_FILES:
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```python {tag}", "```python other"), tag)
        doubled = text + "\n" + f"```python {tag}\n" + extract_tagged_block(text, tag) + "```\n"
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag)
    _p0("[WF-17 s8] (c-page) removed / doubled tags raise LookupError for both blocks")
    for token in ("Study(", ".json", ".npz", ".s4p", ".xdmf"):
        assert token not in text.replace("run_study_file(", ""), token


def test_wf17_s8_malformed_files_refuse(capsys):
    comm = MPI.COMM_WORLD
    check = extract_tagged_block(_text(), CHECK_TAG, lang="bash")
    for name, field in MALFORMED_SAR.items():
        assert f"femem check {name}" in check, name
        got = _field(lambda: load_study(WF1 / name))
        rc = cli.main(["check", str(WF1 / name)])
        rcs = comm.allgather(rc)
        _p0(f"[WF-17 s8] (M) {name}: field {got!r} (want {field!r}); rc per rank {rcs}")
        assert got == field
        assert all(r == 2 for r in rcs)


def test_wf17_s8_records_on_the_page():
    block = extract_tagged_block(_text(), RECORDS_TAG, lang="text")
    got = _records(block)
    _p0(f"[WF-17 s8] (E) records block {got} (item 344: pairs {PAIRS_RECORD}; B1 {B1_RECORD})")
    assert _records_defect(block) is None, _records_defect(block)
    # (c-page), by construction: one digit altered fails (E) and is named
    altered = block.replace(PAIRS_RECORD[1], PAIRS_RECORD[1][:-1] + "5", 1)
    reason = _records_defect(altered)
    _p0(f"[WF-17 s8] (c-page) one digit altered -> {reason}")
    assert reason is not None and PAIRS_LABEL in reason


def test_wf17_s8_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s8] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve windows (-n 2)


def test_wf17_s8_sar_b1_block_reproduces_item344(monkeypatch):
    comm = _solve_preamble()
    t0 = time.perf_counter()
    ns, calls, _ = exec_four_port_block(SAR_TAG, monkeypatch)
    t_block = time.perf_counter() - t0
    run = _check_r_run(ns, calls, SAR_TAG)
    study = run.study
    res = study.fields_result
    iface = study.geometry.sheet_iface
    # the gate's pairing: by solved-sheet azimuth (test_wf1_yaml_sar_b1.py:133-138)
    order = [f"P{s['tag'] - iface}" for s in res.sheets]
    az = {f"P{s['tag'] - iface}": float(s["azimuth_deg"]) for s in res.sheets}
    idx = {pid: _port_index(az[pid], az["P1"], QUADRATURE_STEP_DEG) for pid in az}
    by_k = {k: pid for pid, k in idx.items()}
    assert sorted(by_k) == [0, 1, 2, 3]
    case = {"by_k": by_k, "order": order, "summary": run.summary}

    rel = _s_tie(run)

    # (ii) the 10 g C4 identity, read from the summary
    d = {k: _diag(case, k) for k in range(4)}
    pairs = [abs(d[(k + 1) % 4] - d[k]) / d[k] for k in range(4)]
    pairs_s = [f"{100 * p:.4f}" for p in pairs]
    _p0(
        f"[WF-17 s8] (ii) sheet order {order}; by_k {by_k}; cyclic pairs (%) {pairs_s} "
        f"(record {PAIRS_RECORD}; ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%); {SAR} key"
    )

    # (c1) the far-side ball misses each diagonal
    miss = [abs(_diag(case, k, centre_of=lambda j: (j + 2) % 4) - d[k]) / d[k] for k in range(4)]
    _p0(
        f"[WF-17 s8] (c1) far-side ball miss (%) {[f'{100 * m:.4f}' for m in miss]} "
        f"(RG-SB 88.1162 / 87.9366 / 88.6600 / 88.9315; ASSERTED > "
        f"{C4_COVARIANCE_BAND * 100:.1f}%)"
    )

    # (B1) |B1+| C4 covariance on the run's Study (the gate's images)
    delta = np.radians((az["P2"] - az["P1"]) % 360.0)
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})
    images = {
        "P1@0": ("P1", points),
        "P2@+D": ("P2", _rotate_z(points, delta)),
        "P4@-D": ("P4", _rotate_z(points, -delta)),
        "P3@2D": ("P3", _rotate_z(points, 2.0 * delta)),
    }
    read = {
        label: study.b1_plus_at(pts, drives=[pid], comm=comm)[pid]
        for label, (pid, pts) in images.items()
    }
    n = points.shape[0]
    counts = {label: int(v.sum()) for label, (_, v) in read.items()}
    mask = np.logical_and.reduce([v for _, v in read.values()])
    ref = read["P1@0"][0]
    b1 = [_relative_l2(read[label][0], ref, mask) for label in ("P2@+D", "P4@-D", "P3@2D")]
    b1_s = [f"{100 * v:.4f}" for v in b1]
    _p0(
        f"[WF-17 s8] (B1) valid per image (of {n}) {counts}; P2@+D / P4@-D / P3@2D (%) {b1_s} "
        f"(record {B1_RECORD}; ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%); "
        f"Δ {np.degrees(delta):.6f} deg"
    )
    _p0(
        f"[WF-17 s8] W33 block {t_block:.1f} s (summary elapsed {run.summary['elapsed_s']:.1f} s); "
        f"node {time.perf_counter() - t0:.1f} s"
    )

    assert rel <= LEG_D0_REPRODUCTION_BAND
    for k in range(4):
        assert pairs[k] <= C4_COVARIANCE_BAND, (k, pairs[k])
    assert pairs_s == PAIRS_RECORD, (pairs_s, PAIRS_RECORD)
    assert n > 0
    for label, c in counts.items():
        assert c == n, (label, c, n)
    for v in b1:
        assert v <= C4_COVARIANCE_BAND
    assert b1_s == B1_RECORD, (b1_s, B1_RECORD)
    for k in range(4):
        assert miss[k] > C4_COVARIANCE_BAND, (k, miss[k])


def test_wf17_s8_probes_block_reproduces_u70(monkeypatch):
    _solve_preamble()
    t0 = time.perf_counter()
    ns, calls, _ = exec_four_port_block(PROBES_TAG, monkeypatch)
    t_block = time.perf_counter() - t0
    run = _check_r_run(ns, calls, PROBES_TAG)
    pr = run.summary["probes"]
    res = run.study.fields_result
    sigma, rho = float(run.study.material.sigma), float(run.study.density_kg_m3)
    pts = np.asarray(POINTS)
    order = list(pr["values"])
    sar = np.array([pr["values"][p]["sar_w_per_kg"] for p in order])
    eab = np.array([pr["values"][p]["e_abs_v_per_m"] for p in order])
    hex_equal = [[float(x).hex() for x in p] for p in pr["points"]] == [
        [float(x).hex() for x in p] for p in POINTS
    ]

    rel = _s_tie(run)

    # (U′) the summary is point_sar on the run's own fields, same five-point batch
    own = np.array(
        [
            point_sar(
                res.fields[pid].e_real,
                res.fields[pid].e_imag,
                pts,
                sigma=sigma,
                rho=rho,
                comm=res.comm,
            )
            for pid in order
        ]
    )
    exact = bool(np.array_equal(own, sar))
    ident = float(np.max(np.abs(eab**2 * sigma / (2.0 * rho) - sar) / sar))

    # (c-disp) point 0 displaced 1 mm radially
    disp = np.array(
        [
            point_sar(
                res.fields[pid].e_real,
                res.fields[pid].e_imag,
                np.asarray([DISPLACED_0]),
                sigma=sigma,
                rho=rho,
                comm=res.comm,
            )[0]
            for pid in order
        ]
    )
    c_disp = np.abs(disp - sar[:, 0]) / sar[:, 0]
    _p0(
        f"[WF-17 s8] (PR) ports {order}; points hex equal {hex_equal}; sigma {sigma!r}; "
        f"rho {rho!r}"
    )
    for i, pid in enumerate(order):
        _p0(f"[WF-17 s8]   {pid}: SAR [W/kg per 1 V] {[f'{v:.6e}' for v in sar[i]]}")
    _p0(
        f"[WF-17 s8] (U′) summary == point_sar on run fields (same 5-point batch) {exact}; "
        f"max |E|^2 sigma/(2 rho) vs SAR rel {ident:.3e} (ASSERTED <= {EXACT_IDENTITY_RTOL:.0e}; "
        f"U70 read 2.213e-16)"
    )
    _p0(
        f"[WF-17 s8] (c-disp) per port {[f'{v:.4e}' for v in c_disp]} (U70 6.3865e-02 / "
        f"2.1232e-01 / 8.3482e-02 / 5.6184e-02; ASSERTED > {WIDTH_RTOL:.0e})"
    )
    _p0(
        f"[WF-17 s8] W34 block {t_block:.1f} s (summary elapsed {run.summary['elapsed_s']:.1f} s); "
        f"node {time.perf_counter() - t0:.1f} s"
    )

    assert hex_equal
    assert rel <= LEG_D0_REPRODUCTION_BAND
    assert exact
    assert ident <= EXACT_IDENTITY_RTOL
    assert np.all(sar > 0.0)
    assert float(c_disp.min()) > WIDTH_RTOL
