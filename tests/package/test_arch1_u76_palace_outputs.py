"""`ARCH-1` step 1d (§10 U76): ``femem run --backend palace`` refuses, by name,
every output the Palace route has no path for — the `WF-1` YAML route's
``outputs.sar`` / ``b1_plus`` / ``probes``, ``tune``, ``study.sweep_hz`` and a
structural ``thermal`` — with field ``--backend`` and exit code 2, before any
mesh; the backends page states the carried / refused split from the code.

Brief: docs/planning/chunks/ARCH-1.md, "§9 item 590" and the review ruling
2026-10-22 03:00.  No node builds a mesh or runs Palace.

Anchors: (X) rc 2 by field, nothing written, on 15 rows; (Th) the structural
``thermal`` guard; (A) the admitted path reaches the solve once and carries
``s_matrix`` only; (Pg) the page's carried / refused lines equal the code's;
(Bit) U16's (E2) / (Tie2) / (c-order) / (E1) reading lines string-equal to
U16c's (knob ``FEM_EM_U76_BIT=1``).  Controls (rule (e), asserted, by
construction): (c-drop) the allow-list widened to ``probes`` reproduces the
silent drop; (c-pg) a page copy missing ``outputs.probes`` fails (Pg); (c-bit)
one changed digit fails (Bit).
"""

from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml
from mpi4py import MPI

from fem_em_solver import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    SphereInBoxGeometry,
    StudyConfigError,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    cli,
)
from fem_em_solver import study_yaml
from fem_em_solver.study_yaml import _PALACE_ROUTED_OUTPUTS, _load, _run_palace, _run_study_file_palace

from tests.package.test_wf17_first_study_tutorial import extract_tagged_block
from tests.package.test_wf1_study_yaml import WF1, _fresh_dir, _mapping, _rc_all

REPO = Path(__file__).resolve().parents[2]
SCRATCH = REPO / "logs" / "arch1" / "u76"  # gitignored (logs/)
PAGE = REPO / "docs" / "guide" / "backends.md"
LOG_DIR = REPO / "docs" / "testing" / "logs"
U16C_LOG = "20261018T172228Z_ARCH-1-U16c-W1.log"
BIT_PREFIXES = ("[U16c T4 E2]", "[U16c T4 Tie2]", "[U16c T4 c-order]", "[U16c T5 E1]")
EM_OUTPUT_KEYS = {"s_matrix", "z_matrix", "touchstone", "sar", "b1_plus", "probes",
                  "tune", "temperature", "s_matrices"}
REFUSAL_RE = re.compile(r"palace: (?:a )?(\S+) is not supported")

FOUR_PORT = WF1 / "f_small_four_port.yaml"
TRACKED = ("f_small_em_probes.yaml", "f_small_sar_b1.yaml", "f_small_sweep.yaml",
           "f_small_copper.yaml", "f_small_hole_pec_capacitors.yaml",
           "f_small_hole_pec_tune.yaml", "sphere_mat4.yaml",
           "f_human_tuned_copper_thermal.yaml", "thermal_box_steady_m4.yaml",
           "thermal_ball_wall.yaml")


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _field(fn):
    """(field, message) of the StudyConfigError ``fn`` raises, else (None, None)."""
    try:
        fn()
    except StudyConfigError as exc:
        return exc.field, str(exc)
    return None, None


def _all(comm, flag):
    return bool(comm.allreduce(bool(flag), op=MPI.LAND))


def _stub_builds(monkeypatch):
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called on a refused route")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                ThermalBallGeometry, ThermalBoxGeometry):
        monkeypatch.setattr(cls, "build", stub)
    return calls


def _write(name, mapping, comm):
    path = SCRATCH / "yaml" / name
    if comm.rank == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    comm.Barrier()
    return path


def _rows(comm):
    """(label, file, want field, want key named in the message or None)."""
    base = _mapping(FOUR_PORT)
    sar_b1 = _mapping(WF1 / "f_small_sar_b1.yaml")["outputs"]
    probes = _mapping(WF1 / "f_small_em_probes.yaml")["outputs"]["probes"]
    _fresh_dir(SCRATCH / "yaml", comm)
    rows = [
        ("+outputs.sar", _write("u76_sar.yaml", {**base, "outputs": {"sar": sar_b1["sar"]}}, comm),
         "--backend", "outputs.sar"),
        ("+outputs.b1_plus", _write("u76_b1_plus.yaml",
                                    {**base, "outputs": {"b1_plus": sar_b1["b1_plus"]}}, comm),
         "--backend", "outputs.b1_plus"),
        ("+outputs.probes", _write("u76_probes.yaml",
                                   {**base, "outputs": {"probes": probes}}, comm),
         "--backend", "outputs.probes"),
        ("+tune", _write("u76_tune.yaml", {**base, "tune": {}}, comm), "--backend", "tune"),
    ]
    for name in TRACKED:
        key = "study.sweep_hz" if name == "f_small_sweep.yaml" else None
        rows.append((name, WF1 / name, "--backend", key))
    rows.append(("malformed_thermal_on_port.yaml", WF1 / "malformed_thermal_on_port.yaml",
                 "thermal", None))
    return rows


def _page_lines(text):
    block = extract_tagged_block(text, "palace-outputs", lang="text")
    got = {}
    for line in block.splitlines():
        head, _, rest = line.partition(":")
        got[head.strip()] = rest.split()
    return got


def _page_mismatch(text, honoured, returns, refused):
    """Named differences between the page block and the code's sets ([] = equal)."""
    got = _page_lines(text)
    out = []
    for head, want in (("honoured", honoured), ("returns", returns), ("refused", refused)):
        have = set(got.get(head, []))
        for k in sorted(set(want) - have):
            out.append(f"{head}: missing {k}")
        for k in sorted(have - set(want)):
            out.append(f"{head}: extra {k}")
    return out


# ------------------------------------------------------------------ refusals


def test_u76_palace_refusals_by_field(monkeypatch, capsys):
    """(X), (Th): no mesh, no solve, either build, any width."""
    comm = MPI.COMM_WORLD
    calls = _stub_builds(monkeypatch)
    rows = _rows(comm)
    results = []
    for i, (label, f, want_field, want_key) in enumerate(rows):
        out_py = SCRATCH / "out" / f"row{i:02d}_py"
        out_cli = SCRATCH / "out" / f"row{i:02d}_cli"
        _fresh_dir(out_py, comm)
        _fresh_dir(out_cli, comm)
        g, msg = _field(lambda f=f, o=out_py: _run_study_file_palace(f, output_dir=o, comm=comm))
        capsys.readouterr()
        rc = cli.main(["run", str(f), "--output-dir", str(out_cli), "--backend", "palace"])
        err = capsys.readouterr().err
        rc_all = _rc_all(rc, comm)
        unwritten = _all(comm, not out_py.exists() and not out_cli.exists())
        results.append((label, want_field, want_key, g, msg, rc_all, err, unwritten))
    for label, want_field, want_key, g, msg, rc_all, err, unwritten in results:
        _say(capsys, f"[U76 X] {label}: field {g!r} (want {want_field!r}); msg {msg!r}; "
             f"rc {rc_all}; nothing written {unwritten}")
    for label, want_field, want_key, g, msg, rc_all, err, unwritten in results:
        assert g == want_field, (label, g, msg)
        if want_key is not None:
            assert f"palace: {want_key} is not supported" in msg, (label, msg)
        assert rc_all == (2, 2), (label, rc_all)
        if comm.rank == 0:
            assert "femem: error:" in err, (label, err)
            if want_key is not None:
                assert want_key in err, (label, err)
        assert unwritten, label

    study, kind, outputs, cooling = _load(FOUR_PORT)
    out_th = SCRATCH / "out" / "thermal_guard"
    _fresh_dir(out_th, comm)
    g_th, msg_th = _field(lambda: _run_palace(study, kind, outputs, FOUR_PORT, out_th, "",
                                              0.0, comm, cooling=object()))
    _say(capsys, f"[U76 Th] cooling=object(): field {g_th!r}; msg {msg_th!r}; "
         f"builds {calls['n']}")
    assert g_th == "--backend"
    assert msg_th is not None and msg_th.endswith("palace: thermal is not supported"), msg_th
    assert _all(comm, not out_th.exists())

    builds = comm.allreduce(calls["n"], op=MPI.SUM)
    _say(capsys, f"[U76 X] rows {len(results)}; builds {builds}")
    assert len(results) == 15
    assert builds == 0


# ------------------------------------------------------------------ admitted


def _refused_names(monkeypatch, comm):
    """Every refusal key this module's (X) / (Th) rows raise (rank-local, same on all)."""
    names = set()
    for _label, f, want_field, _k in _rows(comm):
        if want_field != "--backend":
            continue
        _g, msg = _field(lambda f=f: _run_study_file_palace(
            f, output_dir=SCRATCH / "out" / "names", comm=comm))
        names.update(REFUSAL_RE.findall(msg or ""))
    study, kind, outputs, _c = _load(FOUR_PORT)
    _g, msg = _field(lambda: _run_palace(study, kind, outputs, FOUR_PORT,
                                         SCRATCH / "out" / "names", "", 0.0, comm,
                                         cooling=object()))
    names.update(REFUSAL_RE.findall(msg or ""))
    return names


def test_u76_palace_admits_and_carries(monkeypatch, capsys):
    """(A), (Pg), (c-drop), (c-pg): serial by design (``_run_palace`` writes on every rank)."""
    comm = MPI.COMM_WORLD
    if comm.size != 1:
        pytest.skip("serial node: the Palace route refuses a parallel comm by design")
    calls = _stub_builds(monkeypatch)
    spy = {"n": 0}

    def fake_solve(self, comm=MPI.COMM_WORLD, *, backend="dolfinx", **kw):
        if backend != "palace":
            raise AssertionError(f"spy: backend {backend!r}")
        spy["n"] += 1
        return SimpleNamespace(s_matrix=np.eye(4), cells=0, order=1,
                               case_dir=kw.get("case_dir"), sha256={}, palace_total_s=None)

    monkeypatch.setattr("fem_em_solver.study.Study.solve", fake_solve)
    base = _mapping(FOUR_PORT)
    files = [
        ("no outputs", _write("u76_admit_plain.yaml", base, comm)),
        ("outputs.directory", _write(
            "u76_admit_dir.yaml",
            {**base, "outputs": {"directory": str(SCRATCH / "out" / "admit_dir_yaml")}}, comm)),
    ]
    keysets = []
    for label, f in files:
        out = SCRATCH / "out" / f"admit_{f.stem}"
        _fresh_dir(out, comm)
        before = spy["n"]
        run = _run_study_file_palace(f, output_dir=out, comm=comm)
        summary = json.loads((out / "summary.json").read_text())
        keys = set(summary) & EM_OUTPUT_KEYS
        keysets.append(keys)
        _say(capsys, f"[U76 A] {label}: solves {spy['n'] - before}; EM keys {sorted(keys)}; "
             f"run.summary == file {run.summary == summary}")
        assert spy["n"] - before == 1
        assert keys == {"s_matrix"}

    # (c-drop): widen the allow-list to probes -> today's silent drop reproduced
    monkeypatch.setattr(study_yaml, "_PALACE_ROUTED_OUTPUTS", ("directory", "probes"))
    out = SCRATCH / "out" / "c_drop"
    _fresh_dir(out, comm)
    before = spy["n"]
    _run_study_file_palace(WF1 / "f_small_em_probes.yaml", output_dir=out, comm=comm)
    summary = json.loads((out / "summary.json").read_text())
    _say(capsys, f"[U76 c-drop] allow-list + probes: solves {spy['n'] - before}; "
         f"'probes' in summary {'probes' in summary}")
    assert spy["n"] - before == 1
    assert "probes" not in summary
    monkeypatch.setattr(study_yaml, "_PALACE_ROUTED_OUTPUTS", _PALACE_ROUTED_OUTPUTS)

    refused = _refused_names(monkeypatch, comm)
    honoured = ["outputs." + k for k in _PALACE_ROUTED_OUTPUTS]
    returns = sorted(keysets[0])
    text = PAGE.read_text()
    mism = _page_mismatch(text, honoured, returns, refused)
    _say(capsys, f"[U76 Pg] honoured {honoured}; returns {returns}; refused {sorted(refused)}; "
         f"page block {_page_lines(text)}; mismatch {mism}")
    assert mism == [], mism

    # (c-pg): the page with outputs.probes deleted from refused: fails, naming it
    bad = text.replace("refused: outputs.b1_plus outputs.probes ", "refused: outputs.b1_plus ")
    assert bad != text
    mism_c = _page_mismatch(bad, honoured, returns, refused)
    _say(capsys, f"[U76 c-pg] mismatch {mism_c}")
    assert mism_c == ["refused: missing outputs.probes"], mism_c
    assert calls["n"] == 0


# ------------------------------------------------------------------ bitwise


def _bit_knob():
    v = os.environ.get("FEM_EM_U76_BIT")
    if v is None:
        pytest.skip("FEM_EM_U76_BIT unset (W-bit window only)")
    if v != "1":
        raise ValueError(f"FEM_EM_U76_BIT={v!r}: only '1' or unset")


def _read_log(stem):
    plain = LOG_DIR / stem
    if plain.is_file():
        return plain.read_text(errors="replace")
    gz = LOG_DIR / (stem + ".gz")
    if gz.is_file():
        return gzip.decompress(gz.read_bytes()).decode(errors="replace")
    raise FileNotFoundError(stem)


def _bit_lines(text):
    out = {}
    lines = text.splitlines()
    for p in BIT_PREFIXES:
        hits = [ln for ln in lines if ln.startswith(p)]
        if len(hits) != 1:
            raise AssertionError(f"{p}: {len(hits)} lines, need exactly 1")
        out[p] = hits[0]
    return out


def _bit_unequal(a, b):
    return [p for p in BIT_PREFIXES if a[p] != b[p]]


def test_u76_u16_readings_bitwise(capsys):
    """(Bit), (c-bit): log parsing only."""
    _bit_knob()
    ref = _bit_lines(_read_log(U16C_LOG))
    # (c-bit): one changed digit in Tie2 must compare unequal
    tampered = dict(ref)
    tampered["[U16c T4 Tie2]"] = ref["[U16c T4 Tie2]"].replace("2.464e-10", "2.465e-10", 1)
    assert tampered["[U16c T4 Tie2]"] != ref["[U16c T4 Tie2]"]
    c_bit = _bit_unequal(ref, tampered)
    _say(capsys, f"[U76 c-bit] unequal {c_bit}")
    assert c_bit == ["[U16c T4 Tie2]"]

    names = sorted({p.name.removesuffix(".gz") for p in LOG_DIR.iterdir()
                    if re.search(r"_ARCH-1-U76-RGU16[^/]*\.log(\.gz)?$", p.name)})
    assert names, "no ARCH-1-U76-RGU16 log"
    new = _bit_lines(_read_log(names[-1]))
    diff = _bit_unequal(ref, new)
    for p in BIT_PREFIXES:
        _say(capsys, f"[U76 Bit] {p} equal {ref[p] == new[p]}\n  U16c: {ref[p]}\n  U76:  {new[p]}")
    _say(capsys, f"[U76 Bit] comparand {U16C_LOG}; new {names[-1]}; unequal {diff}")
    assert diff == [], diff
