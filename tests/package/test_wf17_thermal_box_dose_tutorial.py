"""`WF-17` step 35 (§9 item 888, §10 U140): the thermal-box page's section
"Read the thermal dose at two points" as a test.

Brief: docs/planning/chunks/WF-17.md "§9 item 888".  The solve node runs in the
**real** build at ``-n 2``, only with ``FEM_EM_WF17_S35=1``; the other nodes
build no mesh.

Anchors (bands imported from item 878's module, none set here): the
``thermal-box-dose`` block, exec'd from ``logs/wf17/s35/``, gives on the R1 file
(I-d) ``cem43_s`` at P5 and ``ADJ_PT`` equal to ``_dose`` of the direct R1
history (``_direct_r1_history``, same window, same two-point batch) at
``IDENTITY_BAND``; (Q-wf) / (Win-wf) through ``_dose_gates`` at ``DOSE_BAR`` /
``AVG_BAR`` on both runs' histories; (E) ≤ ``E_BAND``; (Cells) ``CELLS``.  The
``thermal-box-dose-check`` block's file returns rc 2 and names
``geometry.wall.thickness_m``; (D) 0 dead docrefs, ``MIN_GUIDE_PAGES`` holds.
Controls (asserted, rule (e)): (c-init-d) the R0 run's ``cem43_s`` at P5 misses
R1's by more than ``DOSE_BAR`` (item 878's form and band); (c-page) a removed or
doubled tag raises ``LookupError``.

Not claimed: a dose limit or threshold, a dose under duty or beside probes, a
cited PMMA ρc, the D = 0.25 runs.
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver
from fem_em_solver import cli, load_study
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_yaml_thermal_box_dose import (
    CELLS,
    DOSE_BAR,
    E_BAND,
    HIST_LABELS,
    IDENTITY_BAND,
    MALFORMED,
    R0_DOSE,
    R1_DOSE,
    T_INF,
    _direct_r1_history,
    _dose,
    _dose_gates,
    _field,
)
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf17_thermal_box_wall_tutorial import (
    MIN_GUIDE_PAGES,
    PAGE,
    PAGE_ID,
    WF1_DIR,
    _files_of,
)
from tests.unit.test_doc_reference_exit_codes import run_checker

REPO = PAGE.parents[2]
SCRATCH = REPO / "logs" / "wf17" / "s35"
S35_ENV = "FEM_EM_WF17_S35"
TAG = "thermal-box-dose"
CHECK_TAG = "thermal-box-dose-check"
WANT_FILES = [R1_DOSE.name, R0_DOSE.name]


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S35_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S35_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _text() -> str:
    return PAGE.read_text()


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s35_block_tags_and_files():
    text = _text()
    block = extract_tagged_block(text, TAG)
    files = _files_of(block)
    _p0(
        f"[WF-17 s35] (X) {TAG}: run_study_file( count {block.count('run_study_file(')}; "
        f"FILES {files}; yaml fences {text.count('```yaml')}"
    )
    assert block.count("run_study_file(") == 1 and "Study(" not in block.replace(
        "run_study_file(", ""
    )
    assert files == WANT_FILES
    for name in files:
        assert (WF1_DIR / name).is_file(), name
    assert "```yaml" not in text
    # the s7 ban list, page-wide (`test_wf17_thermal_box_tutorial.py:243–244`)
    bans = {
        name: name in text.replace("run_study_file(", "")
        for name in ("Study(", ".json", ".xdmf", ".npz")
    }
    _p0(f"[WF-17 s35] (X) page-wide bans present {bans}")
    assert not any(bans.values()), bans
    check = extract_tagged_block(text, CHECK_TAG, lang="bash")
    assert check.strip() == f"femem check {MALFORMED.name}"
    # (c-page), asserted: a removed / doubled tag raises
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {TAG}", "```python other"), TAG)
    doubled = text + "\n" + f"```python {TAG}\n" + block + "```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    _p0("[WF-17 s35] (c-page) removed and doubled tag both raise LookupError")


def test_wf17_s35_check_block_refuses():
    comm = MPI.COMM_WORLD
    got = _field(lambda: load_study(MALFORMED))
    rc = cli.main(["check", str(MALFORMED)])
    rcs = comm.allgather(rc)
    _p0(f"[WF-17 s35] (Ref) {MALFORMED.name}: field {got!r}; rc per rank {rcs}")
    assert got == "geometry.wall.thickness_m"
    assert all(r == 2 for r in rcs)


def test_wf17_s35_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s35] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} rc={rc} "
        f"(MIN_GUIDE_PAGES {MIN_GUIDE_PAGES})"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------------ W1 (solve)


def test_wf17_s35_block_reproduces_p68_dose(monkeypatch):
    """W1, real build (the walled-box records are real-build), ``-n 2``: the block's
    R1 dose against the direct R1 history; (Q-wf) / (Win-wf) on both runs."""
    if not _knob():
        pytest.skip(f"{S35_ENV} unset: the walled-box dose block solve is environment-selected")
    comm = MPI.COMM_WORLD
    assert not IS_COMPLEX, "the walled-box records are real-build"
    assert comm.size == 2, comm.size
    t0 = time.perf_counter()

    h_d, _t_d, direct, cells_d, ok_d, s_solve, s_read = _direct_r1_history(comm)
    d_direct = [_dose(h_d[:, j] - T_INF) for j in range(len(HIST_LABELS))]
    _p0(
        f"[WF-17 s35] direct R1: cells {cells_d}; (E) {direct.max_balance_residual:.3e}; "
        f"D {[f'{x:.9e}' for x in d_direct]} s; solve {s_solve:.1f} s, reads {s_read:.1f} s"
    )

    # exec the page block from the scratch (item 838's form)
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), TAG)
    files = _files_of(block)
    if comm.rank == 0:
        if SCRATCH.exists():
            shutil.rmtree(SCRATCH)
        SCRATCH.mkdir(parents=True)
        for name in files:
            (SCRATCH / name).write_bytes((WF1_DIR / name).read_bytes())
    comm.Barrier()
    monkeypatch.chdir(SCRATCH)
    calls = {"n": 0}
    original = fem_em_solver.run_study_file

    def spy(path, *args, **kwargs):
        calls["n"] += 1
        return original(path, *args, **kwargs)

    monkeypatch.setattr(fem_em_solver, "run_study_file", spy)
    ns_ = {"__name__": "__main__"}
    buf = io.StringIO()
    tb = time.perf_counter()
    with contextlib.redirect_stdout(buf):
        exec(compile(block, f"{PAGE_ID}#{TAG}", "exec"), ns_)
    _p0(
        f"[WF-17 s35] block stdout (rank 0): {buf.getvalue()!r}; spy calls {calls['n']}; "
        f"{time.perf_counter() - tb:.1f} s"
    )
    runs = ns_["runs"]
    dz = {"R1": runs[R1_DOSE.name].summary["dose"], "R0": runs[R0_DOSE.name].summary["dose"]}

    gates = {"valid_points": all(ok_d), "runs": len(runs) == 2 and calls["n"] == 2}
    # (I-d)
    idd = [
        abs(dz["R1"]["cem43_s"][j] - d_direct[j]) / abs(d_direct[j])
        for j in range(len(HIST_LABELS))
    ]
    _p0(
        f"[WF-17 s35] (I-d) block R1 cem43_s {[f'{x:.9e}' for x in dz['R1']['cem43_s']]} vs "
        f"direct, rel per point {dict(zip(HIST_LABELS, [f'{x:.3e}' for x in idd]))} "
        f"(IDENTITY_BAND {IDENTITY_BAND:.0e})"
    )
    gates["I_d"] = all(x <= IDENTITY_BAND for x in idd)
    # (Q-wf) / (Win-wf)
    for k in ("R1", "R0"):
        g, rows = _dose_gates(
            k,
            np.asarray(dz[k]["temperature_c"], dtype=float),
            np.asarray(dz[k]["times_s"], dtype=float),
        )
        gates.update(g)
        for r in rows:
            _p0(r)
    # (E) / (Cells)
    e = {k: v.summary["temperature"]["max_balance_residual"] for k, v in runs.items()}
    _p0(
        f"[WF-17 s35] (E) {e} direct {direct.max_balance_residual:.3e} (E_BAND {E_BAND:.0e}); "
        f"(Cells) direct {cells_d} runs {[v.summary['cells'] for v in runs.values()]} "
        f"(CELLS {CELLS})"
    )
    gates["E_direct"] = float(direct.max_balance_residual) <= E_BAND
    gates["E_runs"] = all(v <= E_BAND for v in e.values())
    gates["Cells"] = cells_d == CELLS and all(v.summary["cells"] == CELLS for v in runs.values())
    # (c-init-d), rule (e): asserted at the gate's band, as item 878
    j5 = HIST_LABELS.index("P5")
    d1, d0 = dz["R1"]["cem43_s"][j5], dz["R0"]["cem43_s"][j5]
    miss = abs(d0 / d1 - 1.0)
    _p0(
        f"[WF-17 s35] (c-init-d) P5: R0 D {d0:.6e} s vs R1 D {d1:.6e} s, miss {miss:.4e} "
        f"(asserted > DOSE_BAR {DOSE_BAR:g}); signed R0/R1 - 1 {d0 / d1 - 1.0:+.4e}"
    )
    gates["c_init_d"] = miss > DOSE_BAR
    _p0(f"[WF-17 s35] gates {gates}; elapsed {time.perf_counter() - t0:.1f} s")
    assert all(gates.values()), gates
