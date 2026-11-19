"""`WF-17` step 5 (§9 item 402, chain U45): "the tuned copper birdcage from a
YAML file" is a test.

The committed page `docs/guide/tuned_copper_birdcage_yaml.md` embeds the
tracked `tests/package/wf1/f_human_tuned_copper.yaml` (tag
``f-human-tuned-copper``) and a Python block (tag ``tuned-copper-yaml``) whose
one ``run_study_file`` call is exec'd verbatim from a gitignored scratch
directory under ``logs/``. Asserted (brief: docs/planning/chunks/WF-17.md,
§9 item 402):

(X)  exactly one block per tag (four tags); one ``run_study_file(``, no ``Study(``;
(Y=) the page's YAML block is the tracked file byte for byte;
(K)  configuration, no solve: the loaded Study equals item 382's (``_spec_diff``);
(K-cli) the CLI block makes the Python block's ``run_study_file`` call;
(E)  the output block states item 382b's record at ``.4e``;
(R)  at ``-n 8``, one solve: node S, ``Y_in``, the 18 ccw lattice values,
     SAR10g(c*) and mean |B1+|_ccw reproduce item 390b's comparands / the H9′
     records within ``REUSE_REPRODUCTION_RTOL``;
(D)  ``docrefs --docs-root docs/guide --check-links`` reads dead == 0.

Negative controls: (c1) removed / doubled tags raise; (c2) one-field YAML
mutations diff exactly ``excitation.esr_ohm`` / ``excitation.c_in_f``;
(N2′) the reproduced node S misses the H17 node record (printed factor);
(c3) step 1's planted dead link (run in W0 by node ID).

Disclosed: the solve spy / stop patches the module attribute
``fem_em_solver.study._study_solve_tuned_quadrature`` (item 383 ruling (3);
390b's fixture), the route ``run_study_file`` → ``Study.solve`` takes.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import json
import os
import shlex
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver
import fem_em_solver.study as study_module
import fem_em_solver.study_yaml as study_yaml_module
from fem_em_solver import cli
from fem_em_solver.study_yaml import load_study
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import _spec_diff
from tests.package.test_wf1_yaml_f_human import FH_FILE
from tests.package.test_wf1_yaml_f_human_solve import N_LATTICE, _summary_s
from tests.package.test_wf16_study_tuned_f_human import H17_NODE_JSON, _c_f, _cx, _study
from tests.package.test_wf16_study_tuned_f_human_b1_sar import (
    H9P_CSTAR_LABEL,
    H9P_MEAN_B1_PLUS_CCW_RECORD_T,
    H9P_N_SAMPLE_POINTS,
    H9P_SAR10G_CSTAR_RECORD_W_PER_KG,
    H9P_SAR10G_TABLE_RECORD,
    RECORD_RANK_WIDTH,
)
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_birdcage_b1_plus_map import _sample_points
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_wf10_f_human_tuned_b1_sar import NODE_JSON
from tests.validation.test_wf10_f_human_tuned_sar10g import AZ17_RECORD_DEG, _lattice

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "tuned_copper_birdcage_yaml.md"
SCRATCH = REPO / "logs" / "wf17" / "s5"  # gitignored (/logs/), rank-shared
YAML_NAME = "f_human_tuned_copper.yaml"
YAML_TAG = "f-human-tuned-copper"
PY_TAG = "tuned-copper-yaml"
CLI_TAG = "tuned-copper-yaml-cli"
OUTPUT_TAG = "tuned-copper-yaml-output"
#: Pages in docs/guide's link pass at step 4 (48) plus this one.
MIN_GUIDE_PAGES = 49


class _StopBeforeSolve(Exception):
    pass


class _StopCli(Exception):
    pass


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-17 s5] {msg}", flush=True)


def _page():
    raw = PAGE.read_bytes()
    _p0(f"page sha256 = {_sha(raw)}")
    return raw.decode()


def _blocks(text):
    return {
        "yaml": extract_tagged_block(text, YAML_TAG, lang="yaml"),
        "py": extract_tagged_block(text, PY_TAG),
        "cli": extract_tagged_block(text, CLI_TAG, lang="bash"),
        "out": extract_tagged_block(text, OUTPUT_TAG, lang="text"),
    }


def _call_names(source, name):
    return [
        n
        for n in ast.walk(ast.parse(source))
        if isinstance(n, ast.Call) and getattr(n.func, "id", getattr(n.func, "attr", None)) == name
    ]


def exec_yaml_block(page_text, monkeypatch, *, stop_before_solve):
    """Write the page's YAML into the scratch, chdir there, exec the Python block.

    Returns ``(namespace, recorded, stdout)``. Spies call through and never
    alter what the block computes; with ``stop_before_solve`` the tuned solve
    records its study and raises a sentinel instead.
    """
    comm = MPI.COMM_WORLD
    blocks = _blocks(page_text)
    if comm.rank == 0:
        SCRATCH.mkdir(parents=True, exist_ok=True)
        (SCRATCH / YAML_NAME).write_text(blocks["yaml"])
    comm.Barrier()
    monkeypatch.chdir(SCRATCH)
    recorded = {"solve_calls": 0}
    original_run = fem_em_solver.run_study_file
    real_solve = study_module._study_solve_tuned_quadrature

    def spy_run(path, *args, **kwargs):
        recorded["path"] = path
        recorded["output_dir"] = kwargs.get("output_dir")
        return original_run(path, *args, **kwargs)

    def counting_solve(study, *args, **kwargs):
        recorded["solve_calls"] += 1
        recorded["solve_study"] = study
        if stop_before_solve:
            raise _StopBeforeSolve
        return real_solve(study, *args, **kwargs)

    monkeypatch.setattr(fem_em_solver, "run_study_file", spy_run)
    monkeypatch.setattr(study_module, "_study_solve_tuned_quadrature", counting_solve)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(blocks["py"], f"docs/guide/tuned_copper_birdcage_yaml.md#{PY_TAG}", "exec")
    with contextlib.redirect_stdout(buf):
        try:
            exec(code, namespace)
        except _StopBeforeSolve:
            if not stop_before_solve:
                raise
    out = buf.getvalue()
    _p0(f"block stdout (rank 0): {out!r}")
    return namespace, recorded, out


# ---------------------------------------------------------------- W0 (serial)


def test_wf17_s5_extracts_the_blocks():
    b = _blocks(_page())
    assert all(v.strip() for v in b.values())
    runs = _call_names(b["py"], "run_study_file")
    studies = _call_names(b["py"], "Study")
    _p0(f"(X) run_study_file( calls {len(runs)}, Study( calls {len(studies)}")
    assert len(runs) == 1 and len(studies) == 0


def test_wf17_s5_yaml_block_is_the_tracked_file():
    block = _blocks(_page())["yaml"].encode()
    tracked = FH_FILE.read_bytes()
    _p0(
        f"(Y=) block sha256 {_sha(block)}; tracked sha256 {_sha(tracked)}; "
        f"{len(block)} / {len(tracked)} bytes"
    )
    assert block == tracked


def test_wf17_s5_block_configures_the_h9p_study(monkeypatch):
    text = _page()
    _ns, rec, _out = exec_yaml_block(text, monkeypatch, stop_before_solve=True)
    assert rec["solve_calls"] == 1, "(K) the block never reached the tuned solve"
    resolved = Path(rec["path"]).resolve()
    diff = _spec_diff(rec["solve_study"], _study())
    _p0(
        f"(K) path {rec['path']!r} -> {resolved}; output_dir {rec['output_dir']!r}; "
        f"_spec_diff {diff}; degree {rec['solve_study'].degree}"
    )
    assert resolved == (SCRATCH / YAML_NAME).resolve()
    assert rec["output_dir"] == "femem_out"
    assert diff == []
    assert rec["solve_study"].degree == 1
    assert not (SCRATCH / "femem_out").exists(), "(K) output written before the solve"


def test_wf17_s5_cli_block_matches_the_python_call(monkeypatch):
    b = _blocks(_page())
    line = [ln for ln in b["cli"].splitlines() if "mpiexec -n 8" in ln]
    assert len(line) == 1, line
    tokens = shlex.split(line[0])
    argv = tokens[tokens.index("femem") + 1 :]
    py_call = _call_names(b["py"], "run_study_file")[0]
    py_args = (
        ast.literal_eval(py_call.args[0]),
        ast.literal_eval({k.arg: k.value for k in py_call.keywords}["output_dir"]),
    )
    seen = {}

    def recorder(path, *, output_dir=None, comm=None):
        seen["call"] = (str(path), output_dir)
        raise _StopCli

    monkeypatch.setattr(study_yaml_module, "run_study_file", recorder)
    with pytest.raises(_StopCli):
        cli.main(argv)
    _p0(f"(K-cli) argv {argv}; cli call {seen.get('call')}; python call {py_args}")
    assert seen["call"] == py_args


def test_wf17_s5_planted_page_and_yaml_controls(tmp_path):
    text = _page()
    # (c1) each tag removed, and the YAML block doubled, raise.
    for tag, lang in (
        (YAML_TAG, "yaml"),
        (PY_TAG, "python"),
        (CLI_TAG, "bash"),
        (OUTPUT_TAG, "text"),
    ):
        fence = f"```{lang} {tag}\n"
        assert text.count(fence) == 1, tag
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(fence, f"```{lang}\n"), tag, lang=lang)
    yaml_block = _blocks(text)["yaml"]
    with pytest.raises(LookupError):
        extract_tagged_block(
            text + f"\n```yaml {YAML_TAG}\n{yaml_block}```\n", YAML_TAG, lang="yaml"
        )
    # (c2) one-field mutations: (Y=) fails and _spec_diff names exactly that field.
    tracked = FH_FILE.read_bytes()
    for old, new, field in (
        ("esr_ohm: 0.2", "esr_ohm: 0.0", "excitation.esr_ohm"),
        ("c_in_f: 1.1824370000000001e-10", "c_in_f: 1.182e-10", "excitation.c_in_f"),
    ):
        assert yaml_block.count(old) == 1, old
        mutated = yaml_block.replace(old, new)
        assert mutated.encode() != tracked
        f = tmp_path / f"{field}.yaml"
        f.write_text(mutated)
        diff = _spec_diff(load_study(f), _study())
        _p0(f"(c2) {old!r} -> {new!r}: _spec_diff {diff}")
        assert diff == [field]


def test_wf17_s5_expected_output_states_the_record():
    out = _blocks(_page())["out"]
    expected = f"SAR10g(c*) = {H9P_SAR10G_CSTAR_RECORD_W_PER_KG:.4e} W/kg\n"
    _p0(f"(E) page output {out!r}; expected {expected!r}")
    assert out == expected


def test_wf17_s5_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s5] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- W1 (-n 8)


@complex_only
def test_wf17_s5_block_reproduces_the_u43_record(monkeypatch):
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, f"records are -n {RECORD_RANK_WIDTH}; got {comm.size}"
    # Copied (provenance: tests/package/test_wf1_yaml_f_human_solve.py:127–143,
    # item 390b's ``solved`` fixture): the node JSON read on rank 0 + bcast.
    rec, h17 = None, None
    if comm.rank == 0:
        assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: stop (BLOCKED, rule (d))"
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
        if os.path.exists(H17_NODE_JSON):
            with open(H17_NODE_JSON) as fh:
                h17 = json.load(fh)
    rec, h17 = comm.bcast((rec, h17), root=0)
    assert rec["c_f"] == _c_f(), (rec["c_f"], _c_f())

    text = _page()
    comm.Barrier()
    t0 = time.perf_counter()
    ns, spied, out = exec_yaml_block(text, monkeypatch, stop_before_solve=False)
    comm.Barrier()
    seconds = time.perf_counter() - t0
    run = ns["run"]
    summary = run.summary
    _p0(
        f"block exec (build + node pair + 38 balls + 2 points): {seconds:.1f} s at -n {comm.size}; "
        f"cells {summary.get('cells')}; solve calls {spied['solve_calls']}"
    )

    # (R-1)
    sha_file = _sha(FH_FILE.read_bytes())
    _p0(f"(R-1) summary sha256 {summary['sha256']} vs tracked {sha_file}")
    # (R-S), (R-Y)
    s = _summary_s(summary)
    s_json = np.array([[_cx(v) for v in row] for row in rec["s_node"]], dtype=np.complex128)
    r_s = _residual(s, s_json)
    _p0(
        f"(R-S) s_node vs node JSON residual {r_s:.6e} (band {REUSE_REPRODUCTION_RTOL:g}; "
        f"390b 4.287849e-12)"
    )
    r_y = {}
    for sense in ("ccw", "cw"):
        y = np.array([_cx(p) for p in summary["y_in"][sense]])
        yj = np.array([_cx(p) for p in rec["y_in"][sense]])
        r_y[sense] = float(np.linalg.norm(y - yj) / np.linalg.norm(yj))
        _p0(f"(R-Y) [{sense}] rel {r_y[sense]:.3e} (390b ccw 2.217e-12 / cw 3.216e-12)")
    # (N2′)
    r_n2 = None
    if h17 is not None and "s_node" in h17:
        s17 = np.array([[_cx(v) for v in row] for row in h17["s_node"]], dtype=np.complex128)
        r_n2 = _residual(s, s17)
        _p0(
            f"(N2') s_node vs H17 node at {h17.get('c_f')!r} F residual {r_n2:.6e} = "
            f"{r_n2 / REUSE_REPRODUCTION_RTOL:.3e}x band (> band ASSERTED; prior 4.266875e-02)"
        )
    else:
        _p0("(N2') not run: H17 node JSON absent")
    # (R-SAR)
    labels = [row[0] for row in _lattice(AZ17_RECORD_DEG, 0.0)]
    table = {row[0]: row[1:] for row in H9P_SAR10G_TABLE_RECORD}
    ccw = summary["sar"]["values"]["ccw"]
    cw = summary["sar"]["values"]["cw"]
    sar_c = np.array([d["averaged_sar_w_per_kg"] for d in ccw[:N_LATTICE]])
    r_ccw = np.array([abs(sar_c[i] / table[lab][0] - 1.0) for i, lab in enumerate(labels)])
    key_star = labels[int(np.argmax(sar_c))]
    sar_star = float(ns["sar_star"])
    r_star = abs(sar_star - H9P_SAR10G_CSTAR_RECORD_W_PER_KG) / H9P_SAR10G_CSTAR_RECORD_W_PER_KG
    for i, lab in enumerate(labels):
        _p0(
            f"  [{i:2d}] {lab:<22s} ccw {sar_c[i]:.9e} (rec {table[lab][0]:.9e}, rel {r_ccw[i]:.2e})"
        )
    i_twin = labels.index("r=45 dz=+0 az17+90")
    _p0(
        f"(R-SAR) max rel {r_ccw.max():.3e} (390b 2.857e-10); argmax {key_star} (record "
        f"{H9P_CSTAR_LABEL}); block sar_star {sar_star!r} ({sar_star:.9e}) vs record "
        f"{H9P_SAR10G_CSTAR_RECORD_W_PER_KG:.9e}, rel {r_star:.3e} (390b 4.310e-11)"
    )
    _p0(
        f"printed: cw[18] - cw[{i_twin}] = "
        f"{cw[N_LATTICE]['averaged_sar_w_per_kg'] - cw[i_twin]['averaged_sar_w_per_kg']:.3e}; "
        f"summary b1_plus {summary['b1_plus']['values']}"
    )
    # (R-B)
    res = run.result
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})
    vals, valid = run.study.b1_plus_at(points, drives=("ccw",), comm=comm)["ccw"]
    vals, valid = np.asarray(vals), np.asarray(valid, dtype=bool)
    n_pts, n_valid = int(points.shape[0]), int(np.count_nonzero(valid))
    mean_b = float(np.mean(vals[valid]))
    r_b = abs(mean_b - H9P_MEAN_B1_PLUS_CCW_RECORD_T) / H9P_MEAN_B1_PLUS_CCW_RECORD_T
    _p0(
        f"(R-B) {n_valid}/{n_pts} valid; mean|B1+|_ccw {mean_b:.9e} T vs record "
        f"{H9P_MEAN_B1_PLUS_CCW_RECORD_T:.6e} T, rel {r_b:.3e} (390b 1.101e-07)"
    )
    _p0(f"(R-1) solve calls after every read: {spied['solve_calls']}")

    assert spied["solve_calls"] == 1, spied["solve_calls"]
    assert summary["sha256"] == sha_file
    assert r_s <= REUSE_REPRODUCTION_RTOL, f"(R-S) {r_s:.6e}"
    for sense, r in r_y.items():
        assert r <= REUSE_REPRODUCTION_RTOL, f"(R-Y) {sense} {r:.3e}"
    if r_n2 is not None:
        assert r_n2 > REUSE_REPRODUCTION_RTOL, f"(N2') blind comparator {r_n2:.6e}"
    assert r_ccw.max() <= REUSE_REPRODUCTION_RTOL, f"(R-SAR) {r_ccw.max():.3e}"
    assert key_star == H9P_CSTAR_LABEL, key_star
    assert r_star <= REUSE_REPRODUCTION_RTOL, f"(R-SAR) c* {r_star:.3e}"
    if comm.rank == 0:
        assert f"SAR10g(c*) = {sar_star:.4e} W/kg" in out
    assert n_pts == H9P_N_SAMPLE_POINTS and n_valid == n_pts, (n_pts, n_valid)
    assert r_b <= REUSE_REPRODUCTION_RTOL, f"(R-B) {r_b:.3e}"
