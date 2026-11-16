"""`WF-17` step 25 (§9 item 717, §10 U98): the implant-rod page.

``docs/guide/implant_rod_in_a_study.md`` carries a tagged Python block
(``implant-rod``) and a tagged YAML block (``implant-rod-file``); both run from
the page.

Brief: docs/planning/chunks/WF-17.md, "§9 item 717".

Anchors (imported, unmoved): (Cells) the block's cells == ``ROD_EXACT_CELL_RECORD``;
(S) the block's S vs the direct ``_build_rod`` + ``build_four_port_sweep`` at
``IDENTITY_BAND``, item 709's comparison; (P=Y) Python vs ``femem run`` at
``IDENTITY_BAND``; (D) docrefs 0 dead; (Ref) the page's refusals; (N) the page's
"does not tell you" list.  Printed, never asserted: the tip 1 g Q, a reading with
no record.  Negative control (asserted): the Python block with the rod removed
misses the rod S by more than ``IDENTITY_BAND``; its size is printed against
``POOLED_SEPARATION_FLOOR``.

The solve node runs only with ``FEM_EM_WF17_S25=1`` (then it asserts the complex
build and ``-n 2``, never skips).  Blocks run in a scratch directory under the
gitignored ``logs/``.

`WF-17` step 28 (§9 item 798, §10 U122): the page's "With capacitor sheets"
section.  Its ``implant-rod-caps`` block reproduces item 802's (Id) on the PEC
fixture -- the block's 3x3 against the direct build (``_build_rod`` +
`PORT-24`'s ``_terminated_three_port``, via ``test_wf15_study_rod_caps``'s
``_direct_three_port``) at ``IDENTITY_BAND`` -- with (Cells) ==
``ROD_EXACT_CELL_RECORD``; its inputs are item 802's (``C_F``, ``ESR_OHM``, the
``P2`` sheet).  Control (asserted, no mesh): the page's ``implant-rod-caps-refused``
block raises the degree-2 ``ValueError``, and the copper-hole rod ``Study`` raises
``_refuse_rod_route``'s ``TypeError`` (the one item 802 asserted), each with the
text the page quotes.  The solve node runs only with ``FEM_EM_WF17_S28=1``.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import json
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pytest
import yaml
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, Study
from fem_em_solver.cli import main
from fem_em_solver.study_yaml import StudyConfigError, load_study

from tests.complex_mode import IS_COMPLEX
from tests.mesh.test_birdcage_implant_rod import ROD, _build_rod
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.mesh.test_step_import_four_port import POOLED_SEPARATION_FLOOR
from tests.package.test_wf15_study_implant_rod import ROD_EXACT_CELL_RECORD
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
from tests.validation.test_imp3_palace_twin import _sweep
from tests.validation.test_imp3_rod_birdcage import TIPS
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "implant_rod_in_a_study.md"
PAGE_ID = "docs/guide/implant_rod_in_a_study.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s25"
S25_ENV = "FEM_EM_WF17_S25"
MIN_GUIDE_PAGES = 66

TAG = "implant-rod"
FILE_TAG = "implant-rod-file"
EXPECTED_TAG = "implant-rod-expected"
BANDS_TAG = "implant-rod-bands"
ROD_CALL = "FourPortBirdcageGeometry(implant_rod=rod)"
OFF_CALL = "FourPortBirdcageGeometry()"
S28_ENV = "FEM_EM_WF17_S28"
CAPS_TAG = "implant-rod-caps"
CAPS_EXPECTED_TAG = "implant-rod-caps-expected"
CAPS_BANDS_TAG = "implant-rod-caps-bands"
CAPS_REFUSED_TAG = "implant-rod-caps-refused"
CAPS_REFUSALS_TAG = "implant-rod-caps-refusals"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _text() -> str:
    return PAGE.read_text()


def _knob() -> bool:
    raw = os.environ.get(S25_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S25_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _page_rod(block: str) -> dict:
    for node in ast.parse(block).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "rod":
            return ast.literal_eval(node.value)
    raise LookupError("no `rod = {...}` in the block")


def _page_tips(block: str) -> list:
    for node in ast.parse(block).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "tips":
            return ast.literal_eval(node.value)
    raise LookupError("no `tips = [...]` in the block")


def _scratch(subdir: str) -> Path:
    comm = MPI.COMM_WORLD
    scratch = SCRATCH_ROOT / subdir
    if comm.rank == 0:
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
    comm.Barrier()
    return scratch


def _exec_block(monkeypatch, block: str, label: str):
    scratch = _scratch(label)
    _p0(f"[WF-17 s25] block {label} sha256 {hashlib.sha256(block.encode()).hexdigest()}")
    monkeypatch.chdir(scratch)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(block, f"{PAGE_ID}#{label}", "exec"), namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s25] block {label} stdout (rank 0): {out!r}")
    return namespace, out


def _load_mapping(tmp_path, name, mapping):
    path = tmp_path / f"{name}-{MPI.COMM_WORLD.rank}.yaml"
    path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    return load_study(path)


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s25_blocks_and_refusals(tmp_path):
    from fem_em_solver.backends.palace import solve_palace

    text = _text()
    block = extract_tagged_block(text, TAG)
    mapping = yaml.safe_load(extract_tagged_block(text, FILE_TAG, lang="yaml"))
    expected = extract_tagged_block(text, EXPECTED_TAG, lang="text").strip()
    bands = extract_tagged_block(text, BANDS_TAG, lang="text").strip()

    # The page's inputs are item 709's.
    assert block.count(ROD_CALL) == 1 and block.count("Study(") == 1
    assert _page_rod(block) == ROD
    assert [tuple(t) for t in _page_tips(block)] == list(TIPS.values())
    study = _load_mapping(tmp_path, "page", mapping)
    assert study.geometry == FourPortBirdcageGeometry(implant_rod=ROD)
    assert study.degree == 1 and float(study.frequency_hz) == STEP2D_FREQUENCY_HZ
    assert expected == f"cells = {ROD_EXACT_CELL_RECORD}"
    assert bands == f"Python vs YAML vs the direct build: {IDENTITY_BAND:.0e}"
    _p0(f"[WF-17 s25] (In) rod {_page_rod(block)}, tips {_page_tips(block)}; {expected}; {bands}")

    # (Ref) each refusal the page names.
    cases = {
        "human": (dict(mapping, geometry={"kind": "f_human_birdcage", "as_hole": True}), "implant"),
        "clinical": (dict(mapping, geometry={"kind": "f_clinical_birdcage"}), "implant"),
        "loaded": (dict(mapping, load_result="x.npz", tune={}), "implant"),
    }
    for key in ("radius_m", "length_m", "centre_m", "axis"):
        rod = {k: v for k, v in mapping["implant"]["rod"].items() if k != key}
        cases[f"missing {key}"] = (dict(mapping, implant={"rod": rod}), f"implant.rod.{key}")
    for name, (doc, field) in cases.items():
        with pytest.raises(StudyConfigError) as err:
            _load_mapping(tmp_path, name.replace(" ", "_"), doc)
        _p0(f"[WF-17 s25] (Ref) {name}: {err.value.field}: {err.value.reason}")
        assert err.value.field == field, (name, err.value.field)
    with pytest.raises(ValueError, match="implant_rod"):
        solve_palace(study, MPI.COMM_WORLD, case_dir=tmp_path)

    # (N) the page's limits.
    for token in ("No two-port drive", "No 1 g rung", "No published comparand", "No tissue"):
        assert token in text, token
    for token in ("AED", "HFSS"):
        assert token not in text, token

    # (c-page), by construction: tag removed / doubled raise.
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {TAG}\n", "```python other\n"), TAG)
    with pytest.raises(LookupError):
        extract_tagged_block(text + f"\n```python {TAG}\n{block}```\n", TAG)


def test_wf17_s25_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s25] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s25_implant_rod_blocks(monkeypatch):
    if not _knob():
        pytest.skip(f"{S25_ENV} unset: the page's solve window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port Study needs the complex build"
    assert comm.size == 2, comm.size
    text = _text()
    block = extract_tagged_block(text, TAG)

    t0 = time.perf_counter()
    ns, out = _exec_block(monkeypatch, block, "rod")
    t_rod = time.perf_counter() - t0
    result = ns["result"]
    s_py = np.asarray(result.s_matrix)

    t0 = time.perf_counter()
    scratch = _scratch("yaml")
    file_path = scratch / "implant_rod.yaml"
    if comm.rank == 0:
        file_path.write_text(extract_tagged_block(text, FILE_TAG, lang="yaml"))
    comm.Barrier()
    out_dir = scratch / "out"
    rc = main(["run", str(file_path), "--output-dir", str(out_dir)])
    comm.Barrier()
    summaries = sorted(out_dir.glob("*.json"))
    assert len(summaries) == 1, summaries
    summary = json.loads(summaries[0].read_text())
    s_yaml = np.array([[complex(re, im) for re, im in row] for row in summary["s_matrix"]])
    t_yaml = time.perf_counter() - t0

    t0 = time.perf_counter()
    built = _build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
    direct = _sweep(case, keep=True)
    s_direct = np.asarray(direct.s_matrix)
    t_direct = time.perf_counter() - t0

    t0 = time.perf_counter()
    assert block.count(ROD_CALL) == 1
    ns_off, _ = _exec_block(monkeypatch, block.replace(ROD_CALL, OFF_CALL), "off")
    s_off = np.asarray(ns_off["result"].s_matrix)
    t_off = time.perf_counter() - t0

    def _rel(a, b):
        return comm.allreduce(float(np.linalg.norm(a - b) / np.linalg.norm(b)), op=MPI.MAX)

    rel_direct = _rel(s_py, s_direct)
    rel_yaml = _rel(s_yaml, s_py)
    miss = comm.allreduce(float(np.linalg.norm(s_off - s_py) / np.linalg.norm(s_py)), op=MPI.MIN)
    _p0(f"[WF-17 s25] -n {comm.size}: rod block {t_rod:.1f} s, femem run {t_yaml:.1f} s, direct "
        f"{t_direct:.1f} s, off block {t_off:.1f} s")
    _p0(f"[WF-17 s25] (Cells) block {result.cells}, yaml {summary['cells']}, direct {case['cells']}, "
        f"record {ROD_EXACT_CELL_RECORD}")
    _p0(f"[WF-17 s25] (S) ||S_block - S_direct||/||S_direct|| {rel_direct:.3e}; (P=Y) "
        f"||S_yaml - S_block||/||S_block|| {rel_yaml:.3e} (IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED)")
    _p0(f"[WF-17 s25] (c-off) off block cells {ns_off['result'].cells}; miss {miss:.6e}, "
        f"{miss / IDENTITY_BAND:.3e} x IDENTITY_BAND (POOLED_SEPARATION_FLOOR "
        f"{POOLED_SEPARATION_FLOOR:g}, printed)")
    tips = list(TIPS)
    for pid in sorted(ns["sar"]):
        row = []
        for i, tip in enumerate(tips):
            a = float(ns["sar"][pid][i]["averaged_sar_w_per_kg"])
            b = float(ns_off["sar"][pid][i]["averaged_sar_w_per_kg"])
            row.append(f"{tip} {a / b if b != 0.0 else float('nan'):.6f} (rod {a:.6e}, off {b:.6e} W/kg)")
        _p0(f"[WF-17 s25] tip 1 g Q {pid} (a reading, no record): " + "; ".join(row))

    gates = {
        "stdout": extract_tagged_block(text, EXPECTED_TAG, lang="text").strip() in out.splitlines()
        if comm.rank == 0 else True,
        "cells_record": result.cells == ROD_EXACT_CELL_RECORD,
        "cells_direct": result.cells == case["cells"],
        "cells_yaml": summary["cells"] == ROD_EXACT_CELL_RECORD,
        "S_direct": rel_direct <= IDENTITY_BAND,
        "S_yaml": rel_yaml <= IDENTITY_BAND,
        "rc": comm.allgather(rc) == [0] * comm.size,
        "s4p": (out_dir / "implant_rod.s4p").is_file(),
        "off_misses": miss > IDENTITY_BAND,
        "off_geometry": ns_off["study"].geometry == FourPortBirdcageGeometry(),
    }
    red = sorted(comm.allreduce([k for k, v in gates.items() if not v], op=MPI.SUM))
    _p0(f"[WF-17 s25] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red


# ------------------------------------------------- step 28: the capacitor-sheet section


def _block_study(block: str):
    """The block's ``Study`` built (not solved): its source up to the solve line."""
    head = block.split("result = study.solve()")[0].split("study.solve()")[0]
    namespace = {"__name__": "__main__"}
    exec(compile(head, f"{PAGE_ID}#caps-head", "exec"), namespace)
    return namespace["study"]


class _Reached(Exception):
    pass


def test_wf17_s28_caps_blocks_and_refusals(monkeypatch):
    import fem_em_solver.study as study_mod
    from fem_em_solver import CapacitorSheets, LeontovichWall

    from tests.package.test_wf15_study_implant_rod import _phantom
    from tests.validation.test_ans7_exact_strip import C_F
    from tests.validation.test_port24_series_rc_sheet import ESR_OHM, FREQUENCY_HZ, PLANT_INDEX

    text = _text()
    caps = extract_tagged_block(text, CAPS_TAG)
    refused = extract_tagged_block(text, CAPS_REFUSED_TAG)
    expected = extract_tagged_block(text, CAPS_EXPECTED_TAG, lang="text").strip()
    bands = extract_tagged_block(text, CAPS_BANDS_TAG, lang="text").strip()
    quoted = extract_tagged_block(text, CAPS_REFUSALS_TAG, lang="text").strip().splitlines()

    # Never a mesh in this node: any build raises.
    calls = []

    def _build(self, comm=MPI.COMM_WORLD):
        calls.append(self)
        raise _Reached

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _build)

    # The page's inputs are item 802's.
    all_ports = ("P1", "P2", "P3", "P4")
    sheet = all_ports[PLANT_INDEX]
    st = _block_study(caps)
    assert caps.count(ROD_CALL) == 1 and caps.count("Study(") == 1
    assert _page_rod(caps) == ROD
    assert st.geometry == FourPortBirdcageGeometry(implant_rod=ROD) and not st.geometry.as_hole
    assert st.material == _phantom()
    assert st.excitation == CapacitorSheets(C_F, esr_ohm=ESR_OHM, on=(sheet,))
    assert st.degree == 1 and float(st.frequency_hz) == FREQUENCY_HZ
    assert getattr(st, "conductor", None) is None
    ports = [p for p in all_ports if p != sheet]
    assert expected == f"cells = {ROD_EXACT_CELL_RECORD}; ports = {ports}"
    assert bands == f"Python vs the direct build: {IDENTITY_BAND:.0e}"
    # The refused block is the caps block at degree 2, nothing else.
    assert (refused.replace("degree=2,", "degree=1,").split("study.solve()")[0]
            == caps.split("result = study.solve()")[0])

    # Control: the page's degree-2 block raises the quoted ValueError, before any mesh.
    with pytest.raises(ValueError) as err2:
        exec(compile(refused, f"{PAGE_ID}#{CAPS_REFUSED_TAG}", "exec"), {"__name__": "__main__"})
    # Control: the copper-hole rod Study raises the quoted TypeError (item 802's control 2).
    copper = FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION,
                                      implant_rod=ROD)
    st_cu = Study(copper, _phantom(), CapacitorSheets(C_F, esr_ohm=ESR_OHM, on=(sheet,)),
                  frequency_hz=FREQUENCY_HZ, degree=1,
                  conductor=LeontovichWall(study_mod.COPPER_SIGMA))
    with pytest.raises(TypeError) as err_cu:
        st_cu.solve(MPI.COMM_WORLD)
    _p0(f"\n[WF-17 s28] (In) C_F {C_F!r}, esr {ESR_OHM!r}, sheet {sheet}; {expected}; {bands}")
    _p0(f"[WF-17 s28] (Ref) degree 2: {err2.value}")
    _p0(f"[WF-17 s28] (Ref) copper hole: {err_cu.value}")
    assert quoted == [str(err2.value), str(err_cu.value)], quoted
    assert calls == [], calls

    # (c-page), by construction: the tag removed raises.
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {CAPS_TAG}\n", "```python other\n"),
                             CAPS_TAG)


def test_wf17_s28_caps_block(monkeypatch):
    raw = os.environ.get(S28_ENV, "").strip()
    if raw not in ("", "1"):
        raise ValueError(f"{S28_ENV}={raw!r}: only unset/empty or '1' are accepted")
    if raw != "1":
        pytest.skip(f"{S28_ENV} unset: the caps section's solve window is environment-selected")
    from tests.package.test_wf15_study_rod_caps import _direct_three_port

    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port Study needs the complex build"
    assert comm.size == 2, comm.size
    text = _text()
    caps = extract_tagged_block(text, CAPS_TAG)

    t0 = time.perf_counter()
    ns, out = _exec_block(monkeypatch, caps, "caps")
    t_block = time.perf_counter() - t0
    result = ns["result"]
    s_block = np.asarray(result.s_matrix, dtype=np.complex128)

    t0 = time.perf_counter()
    built = _build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    s_direct, direct_cells = _direct_three_port(built, "rod")
    s_direct = np.asarray(s_direct, dtype=np.complex128)
    t_direct = time.perf_counter() - t0

    rel = comm.allreduce(float(np.linalg.norm(s_block - s_direct) / np.linalg.norm(s_direct)),
                         op=MPI.MAX)
    ports = [p.port_id for p in result.port_definitions]
    _p0(f"[WF-17 s28] -n {comm.size}: caps block {t_block:.1f} s, direct {t_direct:.1f} s")
    _p0(f"[WF-17 s28] (Cells) block {result.cells}, direct {direct_cells}, record "
        f"{ROD_EXACT_CELL_RECORD}; ports {ports}")
    _p0(f"[WF-17 s28] (Id) ||S_block - S_direct||/||S_direct|| {rel:.3e} (IDENTITY_BAND "
        f"{IDENTITY_BAND:.0e}, ASSERTED; item 802 read 4.722e-16)")
    expected = extract_tagged_block(text, CAPS_EXPECTED_TAG, lang="text").strip()
    gates = {
        "stdout": expected in out.splitlines() if comm.rank == 0 else True,
        "cells_record": result.cells == ROD_EXACT_CELL_RECORD,
        "cells_direct": result.cells == direct_cells,
        "shape": s_block.shape == (3, 3) == s_direct.shape,
        "Id": rel <= IDENTITY_BAND,
    }
    red = sorted(comm.allreduce([k for k, v in gates.items() if not v], op=MPI.SUM))
    _p0(f"[WF-17 s28] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red
