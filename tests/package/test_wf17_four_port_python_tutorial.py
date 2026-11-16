"""`WF-17` step 17: the four-port outputs page builds the F-small four-port
`Study` in Python and reads ``Study.sar``, ``Study.b1_plus`` and
``Study.solve_sweep``; it reproduces `WF-16` steps 4 / 5 / 2 (the 10 g C4 pairs,
the ``|B₁⁺|`` C4 covariance, the per-frequency `ANS-4` P1 ties) as the YAML
pages' modules record them.

Brief: docs/planning/chunks/WF-17.md, "§9 items 627, 632, 633, 634" (item 633),
over the template "§9 item 621".  The solve nodes run only with
``FEM_EM_WF17_S17`` set to ``sar`` / ``b1`` / ``sweep`` (each node runs on its
own value and skips on another; any other value raises), and then assert the
complex build and ``-n ANS4_EXACT_STRIP_RECORD_RANK_WIDTH``, never skip.  The
other nodes build no mesh.

Anchors: (X) one block per tag (four tags); each python block holds one
``Study(`` call and exactly its own output call(s) by AST name, and no
``run_study_file``; (In) each block's Study equals the gate's inputs exactly
(no mesh; the output calls patched to record); (E) the bands block = the
imported bands at ``.0e``; (R) the S tie at ``LEG_D0_REPRODUCTION_BAND``;
(ii) the 10 g cyclic pairs at ``C4_COVARIANCE_BAND`` and at
``PAIRS_RECORD``'s ``.4f`` digits; (B1) the ``|B₁⁺|`` covariance at
``C4_COVARIANCE_BAND`` and at ``B1_RECORD``'s digits (both records imported by
name from the YAML page's module and compared as it compares them); (R′) the
three sweep ties at ``LEG_D0_REPRODUCTION_BAND`` on one mesh of
``ANS4_EXACT_STRIP_CELL_RECORD`` cells, band only, digits printed beside
``W1A_R`` (the 64 MHz digit differs between routes); (D) 0 dead links.
Controls (rule (e)), asserted and backed: (c1) the far-side ball misses each
diagonal (U73 88.1162 / 87.9366 / 88.6600 / 88.9315 %); (c2) P2 at −Δ misses
(item 322 35.0040 %); (c-swap) the 10 ↔ 128 MHz exchange misses (W1a 1.2367 /
1.2763); (c-page) by construction.

Mesh economy (disclosed, as the template's): inside each solve node
``FourPortBirdcageGeometry.build`` is wrapped so the first call builds the real
mesh and later equal-geometry calls on the same communicator return it.  The
blocks are exec'd from the gitignored scratch ``logs/wf17/s17/``.
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
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf16_study_sar import _phantom
from tests.package.test_wf16_study_sweep import C4_HAND_TABLE, KEYS, _records
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf17_four_port_outputs_tutorial import B1_RECORD, PAIRS_RECORD
from tests.package.test_wf17_sweep_tutorial import W1A_R
from tests.package.test_wf1_yaml_sar_b1 import SAR, _diag
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    FREQUENCIES,
)
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import QUADRATURE_STEP_DEG, _port_index
from tests.validation.test_birdcage_sar_mass_averaged import CENTRE_RADIUS_M
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "four_port_outputs_from_python.md"
PAGE_ID = "docs/guide/four_port_outputs_from_python.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s17"
S17_ENV = "FEM_EM_WF17_S17"
S17_VALUES = ("sar", "b1", "sweep")
MIN_GUIDE_PAGES = 66

SAR_TAG = "sar-python"
B1_TAG = "b1-python"
SWEEP_TAG = "sweep-python"
BANDS_TAG = "four-port-python-bands"
# per block: the output calls it must hold (once each), and the ones it must not
OWN_CALLS = {
    SAR_TAG: ("solve", "sar"),
    B1_TAG: ("b1_plus", "b1_plus_at"),
    SWEEP_TAG: ("solve_sweep",),
}
ALL_CALLS = ("solve", "sar", "b1_plus", "b1_plus_at", "solve_sweep", "run_study_file", "tune")

# Printed beside the readings; never asserted (the backing of each control).
U73_C1 = "88.1162 / 87.9366 / 88.6600 / 88.9315"  # `20261021T003305Z_WF-17-s8-W33.log:3625`
ITEM322_C2 = "35.0040"  # test_wf1_yaml_sar_b1.py's (c2) note
W1A_SWAP = "1.2367 / 1.2763"  # test_wf17_sweep_tutorial.py's (c-swap) note


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


def _block_defect(block: str, tag: str):
    """None when the python block has its (X) shape, else the reason."""
    n_study = len(_calls_named(block, "Study"))
    if n_study != 1:
        return f"{n_study} Study calls; need exactly one"
    for name in ALL_CALLS:
        n = len([c for c in _calls_named(block, name)])
        want = 1 if name in OWN_CALLS[tag] else 0
        if n != want:
            return f"{n} {name} calls; need {want}"
    return None


def _bands_lines() -> list:
    return [
        f"S vs the P1 records, per frequency: {LEG_D0_REPRODUCTION_BAND:.0e}",
        f"10 g C4 cyclic pairs: {C4_COVARIANCE_BAND:.0e}",
        f"|B1+| C4 covariance: {C4_COVARIANCE_BAND:.0e}",
    ]


def _bands_defect(block: str):
    got = block.strip("\n").splitlines()
    want = _bands_lines()
    if len(got) != len(want):
        return f"{len(got)} lines != {len(want)}"
    for g, w in zip(got, want):
        if g != w:
            return f"line {g!r} != {w!r}"
    return None


def _knob():
    raw = os.environ.get(S17_ENV, "").strip()
    if raw == "":
        return None
    if raw in S17_VALUES:
        return raw
    raise ValueError(f"{S17_ENV}={raw!r}: only unset/empty or one of {S17_VALUES} are accepted")


def _solve_preamble(which: str):
    got = _knob()
    if got != which:
        pytest.skip(f"{S17_ENV}={got!r}: this window is {which!r} (complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port port Study needs the complex build"
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def _prepare_scratch():
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    if comm.rank == 0:
        print(f"[WF-17 s17] page sha256 {_sha(raw)}", flush=True)
        if SCRATCH_ROOT.exists():
            shutil.rmtree(SCRATCH_ROOT)
        SCRATCH_ROOT.mkdir(parents=True)
    comm.Barrier()


def exec_block(monkeypatch, tag):
    """chdir into the scratch on every rank, exec the tagged block as __main__."""
    block = extract_tagged_block(PAGE.read_text(), tag)
    _p0(f"[WF-17 s17] block {tag} sha256 {_sha(block.encode())}")
    monkeypatch.chdir(SCRATCH_ROOT)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{tag}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s17] block {tag} stdout (rank 0): {out!r}")
    return namespace, out


def _assert_gate_inputs(st, label):
    _p0(f"[WF-17 s17] (In) {label}: geometry == f_small() {st.geometry == FourPortBirdcageGeometry.f_small()}; "
        f"material {st.material!r}; excitation {st.excitation!r}; frequency "
        f"{float(st.frequency_hz).hex()} (want {FREQUENCY_64_HZ.hex()}); degree {st.degree}; "
        f"conductor {st.conductor!r}")
    assert type(st.geometry) is FourPortBirdcageGeometry
    assert st.geometry == FourPortBirdcageGeometry.f_small()
    assert st.material == _phantom()
    assert float(st.material.sigma).hex() == SALINE_SIGMA.hex()
    assert float(st.material.epsilon_r).hex() == SALINE_EPSILON_R.hex()
    assert type(st.excitation) is LumpedPortDrive and st.excitation == LumpedPortDrive()
    assert float(st.frequency_hz).hex() == FREQUENCY_64_HZ.hex()
    assert st.degree == 1
    assert st.conductor is None


def _run_block(tag, monkeypatch, comm):
    """The template's mesh economy: one real build, equal-geometry calls hit."""
    real_build = FourPortBirdcageGeometry.build
    cache = {"key": None, "built": None, "hits": 0, "real": 0}

    def _build(self, comm_=MPI.COMM_WORLD):
        if cache["built"] is not None and cache["key"] == self and comm_ is comm:
            cache["hits"] += 1
            return cache["built"]
        cache["real"] += 1
        out = real_build(self, comm_)
        cache["key"], cache["built"] = self, out
        return out

    _prepare_scratch()
    t0 = time.perf_counter()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        ns, out = exec_block(monkeypatch, tag)
    t_block = time.perf_counter() - t0
    _p0(f"[WF-17 s17] build cache: {cache['real']} real build, {cache['hits']} hits; block {tag} "
        f"{t_block:.1f} s")
    assert cache["real"] == 1, cache["real"]
    _assert_gate_inputs(ns["study"], tag)
    return ns


def _pairing(study):
    res = study.fields_result
    iface = study.geometry.sheet_iface
    order = [f"P{s['tag'] - iface}" for s in res.sheets]
    az = {f"P{s['tag'] - iface}": float(s["azimuth_deg"]) for s in res.sheets}
    idx = {pid: _port_index(az[pid], az["P1"], QUADRATURE_STEP_DEG) for pid in az}
    by_k = {k: pid for pid, k in idx.items()}
    assert sorted(by_k) == [0, 1, 2, 3]
    return res, order, az, by_k


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s17_blocks_and_inputs(monkeypatch, tmp_path):
    text = _text()
    blocks = {tag: extract_tagged_block(text, tag) for tag in OWN_CALLS}
    extract_tagged_block(text, BANDS_TAG, lang="text")
    for tag, block in blocks.items():
        reason = _block_defect(block, tag)
        counts = {n: len(_calls_named(block, n)) for n in ("Study",) + ALL_CALLS}
        _p0(f"[WF-17 s17] (X) {tag}: calls {counts}; defect {reason}")
        assert reason is None, reason
        assert block.count("Study(") == 1 and "run_study_file" not in block
    for link in ("(four_port_outputs_from_yaml.md)", "(sweep_from_yaml.md)",
                 "(../../tests/package/test_wf17_four_port_python_tutorial.py)"):
        assert link in text, link

    # (c-page), by construction: a removed / doubled tag raises; a block without
    # its output call, or with another block's, fails (X) naming the call
    for tag, lang in [(t, "python") for t in OWN_CALLS] + [(BANDS_TAG, "text")]:
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    planted = tmp_path / "page.md"
    planted.write_text(text.replace("sar = study.sar(", "sar = study.b1_plus(", 1))
    r_sar = _block_defect(extract_tagged_block(planted.read_text(), SAR_TAG), SAR_TAG)
    planted.write_text(text.replace("study.solve_sweep(", "study.solve(", 1))
    r_sweep = _block_defect(extract_tagged_block(planted.read_text(), SWEEP_TAG), SWEEP_TAG)
    _p0(f"[WF-17 s17] (c-page) tags removed / doubled raise LookupError for all four blocks; "
        f"sar -> b1_plus: {r_sar}; solve_sweep -> solve: {r_sweep}")
    assert r_sar is not None and ("sar" in r_sar or "b1_plus" in r_sar)
    assert r_sweep is not None and ("solve" in r_sweep)
    for token in (".json", ".npz", ".s4p", ".xdmf", "AED", "HFSS"):
        assert token not in text, token

    # (In) each block's Study, its output calls patched to record (no mesh)
    calls = []
    fake_sheets = [{"tag": 0, "azimuth_deg": 0.0}]

    def rec(name, value):
        def f(self, *args, **kwargs):
            calls.append((name, self, args, kwargs))
            return value
        return f

    sweep_stub = SimpleNamespace(frequencies_hz=(), s_matrices=np.zeros((0, 4, 4)))
    monkeypatch.setattr(Study, "solve", rec("solve", SimpleNamespace(sheets=fake_sheets)))
    monkeypatch.setattr(Study, "sar", rec("sar", {"P1": [{}]}))
    monkeypatch.setattr(Study, "b1_plus", rec("b1_plus", {"P1": None}))
    monkeypatch.setattr(Study, "b1_plus_at",
                        rec("b1_plus_at", {"P1": (np.zeros(1), np.ones(1, dtype=bool))}))
    monkeypatch.setattr(Study, "solve_sweep", rec("solve_sweep", sweep_stub))
    _prepare_scratch()
    for tag in OWN_CALLS:
        calls.clear()
        ns, _ = exec_block(monkeypatch, tag)
        st = ns["study"]
        names = [c[0] for c in calls]
        _p0(f"[WF-17 s17] (In) {tag}: calls {[(c[0], c[2], c[3]) for c in calls]}")
        assert names == list(OWN_CALLS[tag]), names
        assert all(c[1] is st for c in calls)
        _assert_gate_inputs(st, tag)
        if tag == SAR_TAG:
            sar_call = calls[1]
            assert sar_call[3]["mass_g"] == 10.0 and "radius_m" not in sar_call[3]
            assert "drives" not in sar_call[3]
        if tag == B1_TAG:
            assert calls[0][2] == () and calls[0][3] == {}
            assert [tuple(p) for p in calls[1][2][0]] == [(0.0, 0.0, 0.0)]
        if tag == SWEEP_TAG:
            got = [float(f).hex() for f in calls[0][2][0]]
            want = [float(f).hex() for f in FREQUENCIES.values()]
            _p0(f"[WF-17 s17] (In) sweep frequencies hex {got} (want {want})")
            assert got == want


def test_wf17_s17_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s17] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _bands_defect(block) is None, _bands_defect(block)
    altered = block.replace(f"{LEG_D0_REPRODUCTION_BAND:.0e}",
                            f"{LEG_D0_REPRODUCTION_BAND * 10:.0e}", 1)
    reason = _bands_defect(altered)
    _p0(f"[WF-17 s17] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "P1 records" in reason


def test_wf17_s17_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s17] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve windows (-n 2)


def test_wf17_s17_sar_block_reproduces_wf16_s4(monkeypatch):
    comm = _solve_preamble("sar")
    t0 = time.perf_counter()
    ns = _run_block(SAR_TAG, monkeypatch, comm)
    study, result, sar, centres = ns["study"], ns["result"], ns["sar"], ns["centres"]
    res, order, az, by_k = _pairing(study)

    # the block's centres = the gate's expression on the solved sheets, by float.hex
    own = [[float(x).hex() for x in (CENTRE_RADIUS_M * float(np.cos(np.radians(az[pid]))),
                                     CENTRE_RADIUS_M * float(np.sin(np.radians(az[pid]))), 0.0)]
           for pid in order]
    got = [[float(x).hex() for x in c] for c in centres]
    same_sheets = [s["tag"] for s in result.sheets] == [s["tag"] for s in res.sheets]

    # (R) the block's solve() S vs the degree-1 record
    rel = _rel(np.asarray(result.s_matrix, dtype=np.complex128), _record(1))
    _p0(f"[WF-17 s17] (R) cells {result.cells}; ||S - P1|| / ||P1|| = {rel:.3e} (ASSERTED <= "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}); centres hex = gate expression {got == own}; sheet "
        f"order solve() = fields {same_sheets}")

    # (ii) the 10 g C4 identity, the YAML module's comparison on the block's table
    case = {"by_k": by_k, "order": order, "summary": {"sar": {"values": sar}}}
    d = {k: _diag(case, k) for k in range(4)}
    pairs = [abs(d[(k + 1) % 4] - d[k]) / d[k] for k in range(4)]
    pairs_s = [f"{100 * p:.4f}" for p in pairs]
    _p0(f"[WF-17 s17] (ii) sheet order {order}; by_k {by_k}; cyclic pairs (%) {pairs_s} "
        f"(PAIRS_RECORD {PAIRS_RECORD}; ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}% and digits "
        f"equal); {SAR} key")

    # (c1) the far-side ball misses each diagonal
    miss = [abs(_diag(case, k, centre_of=lambda j: (j + 2) % 4) - d[k]) / d[k] for k in range(4)]
    _p0(f"[WF-17 s17] (c1) far-side ball miss (%) {[f'{100 * m:.4f}' for m in miss]} "
        f"(U73 {U73_C1}; ASSERTED > {C4_COVARIANCE_BAND * 100:.1f}%)")
    _p0(f"[WF-17 s17] W-sar node {time.perf_counter() - t0:.1f} s")

    assert result.cells == ANS4_EXACT_STRIP_CELL_RECORD
    assert got == own
    assert same_sheets
    assert sorted(sar) == ["P1", "P2", "P3", "P4"]
    assert all(len(v) == 4 for v in sar.values())
    assert rel <= LEG_D0_REPRODUCTION_BAND
    for k in range(4):
        assert pairs[k] <= C4_COVARIANCE_BAND, (k, pairs[k])
    assert pairs_s == PAIRS_RECORD, (pairs_s, PAIRS_RECORD)
    for k in range(4):
        assert miss[k] > C4_COVARIANCE_BAND, (k, miss[k])


def test_wf17_s17_b1_block_reproduces_wf16_s5(monkeypatch):
    comm = _solve_preamble("b1")
    t0 = time.perf_counter()
    ns = _run_block(B1_TAG, monkeypatch, comm)
    study, maps, b1 = ns["study"], ns["maps"], ns["b1"]
    res, order, az, by_k = _pairing(study)

    rel = _rel(np.asarray(res.s_matrix, dtype=np.complex128), _record(1))
    _p0(f"[WF-17 s17] (R) cells {res.cells}; ||S_keep - P1|| / ||P1|| = {rel:.3e} (ASSERTED <= "
        f"{LEG_D0_REPRODUCTION_BAND:.0e})")

    # (B1) the YAML module's images on the block's Study
    delta = np.radians((az["P2"] - az["P1"]) % 360.0)
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})
    images = {
        "P1@0": ("P1", points),
        "P2@+D": ("P2", _rotate_z(points, delta)),
        "P4@-D": ("P4", _rotate_z(points, -delta)),
        "P3@2D": ("P3", _rotate_z(points, 2.0 * delta)),
    }
    read = {label: study.b1_plus_at(pts, drives=[pid], comm=comm)[pid]
            for label, (pid, pts) in images.items()}
    control = study.b1_plus_at(_rotate_z(points, -delta), drives=["P2"], comm=comm)["P2"]
    n = points.shape[0]
    counts = {label: int(v.sum()) for label, (_, v) in read.items()}
    mask = np.logical_and.reduce([v for _, v in read.values()] + [control[1]])
    ref = read["P1@0"][0]
    vals = [_relative_l2(read[label][0], ref, mask) for label in ("P2@+D", "P4@-D", "P3@2D")]
    vals_s = [f"{100 * v:.4f}" for v in vals]
    mask4 = np.logical_and.reduce([v for _, v in read.values()])
    v4 = [_relative_l2(read[label][0], ref, mask4) for label in ("P2@+D", "P4@-D", "P3@2D")]
    vals4_s = [f"{100 * v:.4f}" for v in v4]
    miss = _relative_l2(control[0], ref, mask)
    _p0(f"[WF-17 s17] (B1) valid per image (of {n}) {counts}, control {int(control[1].sum())}; "
        f"P2@+D / P4@-D / P3@2D (%) {vals4_s} (B1_RECORD {B1_RECORD}; ASSERTED <= "
        f"{C4_COVARIANCE_BAND * 100:.1f}% and digits equal); with the control in the mask "
        f"{vals_s} (printed); Δ {np.degrees(delta):.6f} deg; block maps {sorted(maps)}, centre "
        f"valid {[bool(v.all()) for _, v in b1.values()]}")
    _p0(f"[WF-17 s17] (c2) P2@-D vs P1: {miss * 100:.4f}% (item 322 {ITEM322_C2}%; ASSERTED > "
        f"{C4_COVARIANCE_BAND * 100:.1f}%)")
    _p0(f"[WF-17 s17] W-b1 node {time.perf_counter() - t0:.1f} s")

    assert res.cells == ANS4_EXACT_STRIP_CELL_RECORD
    assert sorted(maps) == ["P1", "P2", "P3", "P4"]
    assert sorted(b1) == ["P1", "P2", "P3", "P4"]
    assert all(bool(v.all()) for _, v in b1.values())
    assert rel <= LEG_D0_REPRODUCTION_BAND
    assert n > 0
    for label, c in counts.items():
        assert c == n, (label, c, n)
    assert int(control[1].sum()) == n
    for v in v4:
        assert v <= C4_COVARIANCE_BAND
    assert vals4_s == B1_RECORD, (vals4_s, B1_RECORD)
    assert miss > C4_COVARIANCE_BAND


def test_wf17_s17_sweep_block_reproduces_wf16_s2(monkeypatch):
    comm = _solve_preamble("sweep")
    t0 = time.perf_counter()
    ns = _run_block(SWEEP_TAG, monkeypatch, comm)
    study, sw = ns["study"], ns["sweep"]
    rec = _records(1)
    s_all = np.asarray(sw.s_matrices, dtype=np.complex128)
    rels = {k: _rel(s_all[i], rec[k]) for i, k in enumerate(KEYS)}
    one_mesh = all(r.mesh is sw.results[0].mesh for r in sw.results)
    hex_ok = [float(f).hex() for f in sw.frequencies_hz] == [float(f).hex()
                                                            for f in FREQUENCIES.values()]
    i10, i128 = KEYS.index("10"), KEYS.index("128")
    sw_10 = _rel(s_all[i10], rec["128"])
    sw_128 = _rel(s_all[i128], rec["10"])
    _p0(f"[WF-17 s17] (R') ||S_f - P1_f||/||P1_f|| "
        + " / ".join(f"{k} MHz {v:.3e}" for k, v in rels.items())
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}; digits printed, never asserted; "
        f"YAML route W1a {W1A_R}); cells {sw.results[0].cells} (record "
        f"{ANS4_EXACT_STRIP_CELL_RECORD}); one mesh {one_mesh}; frequencies hex equal {hex_ok}; "
        f"study.result is None {study.result is None}")
    _p0(f"[WF-17 s17] (c-swap) S_10 vs P1_128 {sw_10:.4f} (hand {C4_HAND_TABLE[1][('10', '128')]}); "
        f"S_128 vs P1_10 {sw_128:.4f} (hand {C4_HAND_TABLE[1][('128', '10')]}) "
        f"(ASSERTED > {LEG_D0_REPRODUCTION_BAND:.0e}; W1a {W1A_SWAP})")
    _p0(f"[WF-17 s17] W-sweep mesh {sw.mesh_time:.2f} s; sweeps "
        f"{[round(t, 2) for t in sw.sweep_times]}; node {time.perf_counter() - t0:.1f} s")

    assert len(sw.results) == len(KEYS)
    assert s_all.shape == (len(KEYS), 4, 4)
    assert sw.results[0].cells == ANS4_EXACT_STRIP_CELL_RECORD
    assert one_mesh
    assert hex_ok
    assert study.result is None
    for k, v in rels.items():
        assert v <= LEG_D0_REPRODUCTION_BAND, k
    assert sw_10 > LEG_D0_REPRODUCTION_BAND
    assert sw_128 > LEG_D0_REPRODUCTION_BAND
