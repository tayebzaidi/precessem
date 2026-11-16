"""`WF-17` step 14: the capacitor-sheets-with-ESR tutorial page builds the copper
F-small ``Study`` with a 0.2 ohm series-RC sheet on P2 and reproduces `WF-16`
step 10's per-sheet circuit-reduction identity, with the unterminated 4x4's
copper tie to the `ANS-6` degree-1 record.

Brief: docs/planning/chunks/WF-17.md, "§9 item 621".  The solve node runs only
with ``FEM_EM_WF17_S14=1`` (and then asserts the complex build and ``-n
ANS6_EXACT_STRIP_RECORD_RANK_WIDTH``, never skips).  The other nodes build no
mesh.

Anchors: (X) one block per tag; the python block holds one ``Study(`` call, one
``CapacitorSheets(`` call with keywords ``esr_ohm`` and ``on``, one ``.solve(``
call and no ``run_study_file`` call (AST names); (In) the block's Study inputs
equal the gate's exactly (``Study.solve`` patched, no mesh); (E) the bands block
= the three imported constants at ``.0e``; (R0) S0's cells and class means vs
the copper record at ``LEG_D0_REPRODUCTION_BAND``; (I) the page's 3x3 vs the
circuit reduction of S0 with Z = r + 1/(jwC) on P2 at ``REDUCTION_BAND``; (D) 0
dead links.  Controls (rule (e)): (N1) the planted 0.4 ohm on P2 misses the
0.2 ohm reduction (asserted, backed by `WF-16` step 10's 10.953x); (c-page) by
construction.  Bands only are asserted; every reading is printed.

Mesh economy (disclosed, as the gate's ``s10`` fixture): inside the solve node
``FourPortBirdcageGeometry.build`` is wrapped so the first call builds the real
mesh and later equal-geometry calls on the same communicator return it.  The
(In) node's ``Study.solve`` patch returns a stub carrying the two attributes the
block prints (``port_definitions``, ``capacitor_impedance_ohm``) instead of
``None``, so the block's print line runs without a mesh (disclosed adaptation).
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

import pytest
from mpi4py import MPI

from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.ports.circuit import reduce_terminated_ports
from fem_em_solver.study import COPPER_SIGMA, GATED_WIDTH_FRACTION

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_study_capacitor_sheets import PORTS, _phantom, _study
from tests.package.test_wf16_study_copper import _geometry
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_ans7_exact_strip import C_F
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port24_series_rc_sheet import (
    ESR_OHM,
    PLANT_INDEX,
    PLANT_OHM,
    _rc,
    _rel,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "capacitor_sheets_with_esr.md"
PAGE_ID = "docs/guide/capacitor_sheets_with_esr.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s14"
S14_ENV = "FEM_EM_WF17_S14"
MIN_GUIDE_PAGES = 63

ESR_TAG = "capacitor-esr"
BANDS_TAG = "capacitor-esr-bands"

# Printed beside the readings (`20261011T202050Z_WF-16-s10-W1.log` :976 / :980
# / :996); never asserted.
W1_R0 = "4.093e-11 / 8.673e-11 / 1.376e-10"
W1_I_P2 = "3.036163e-05"
W1_N1 = "1.095344e-02 = 10.953x"


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


def _block_defect(block: str):
    """None when the python block has the (X) shape, else the reason."""
    n_study = len(_calls_named(block, "Study"))
    if n_study != 1:
        return f"{n_study} Study calls; need exactly one"
    caps = _calls_named(block, "CapacitorSheets")
    if len(caps) != 1:
        return f"{len(caps)} CapacitorSheets calls; need exactly one"
    keywords = {k.arg for k in caps[0].keywords}
    for kw in ("esr_ohm", "on"):
        if kw not in keywords:
            return f"CapacitorSheets call lacks the keyword {kw!r}"
    solves = [c for c in _calls_named(block, "solve") if isinstance(c.func, ast.Attribute)]
    if len(solves) != 1:
        return f"{len(solves)} .solve( calls; need exactly one"
    if _calls_named(block, "run_study_file"):
        return "the block calls run_study_file"
    return None


def _bands_lines() -> list:
    return [
        f"Study S vs the circuit reduction with Z = r + 1/(jωC): {REDUCTION_BAND:.0e}",
        f"Unterminated 4×4 vs the copper record: {LEG_D0_REPRODUCTION_BAND:.0e}",
        f"Series resistance r (Ω): {ESR_OHM:.0e}",
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


def _knob() -> bool:
    raw = os.environ.get(S14_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S14_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(f"{S14_ENV} unset: the tutorial solve window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the copper four-port capacitor-sheet Study needs the complex build"
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def exec_esr_block(monkeypatch):
    """Clear the scratch, chdir there on every rank, exec the block as __main__."""
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), ESR_TAG)
    if comm.rank == 0:
        print(f"[WF-17 s14] page sha256 {_sha(raw)}; block {ESR_TAG} sha256 "
              f"{_sha(block.encode())}", flush=True)
        if SCRATCH_ROOT.exists():
            shutil.rmtree(SCRATCH_ROOT)
        SCRATCH_ROOT.mkdir(parents=True)
    comm.Barrier()
    monkeypatch.chdir(SCRATCH_ROOT)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{ESR_TAG}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s14] block {ESR_TAG} stdout (rank 0): {out!r}")
    return namespace, out


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s14_blocks_and_inputs(monkeypatch, tmp_path):
    text = _text()
    block = extract_tagged_block(text, ESR_TAG)
    reason = _block_defect(block)
    _p0(f"[WF-17 s14] (X) {ESR_TAG}: Study {len(_calls_named(block, 'Study'))}, "
        f"CapacitorSheets {len(_calls_named(block, 'CapacitorSheets'))}, solve "
        f"{len(_calls_named(block, 'solve'))}, run_study_file "
        f"{len(_calls_named(block, 'run_study_file'))}; defect {reason}")
    assert reason is None, reason
    assert block.count("Study(") == 1 and block.count("run_study_file") == 0
    extract_tagged_block(text, BANDS_TAG, lang="text")

    # (c-page), by construction: a removed / doubled tag raises
    for tag, lang in ((ESR_TAG, "python"), (BANDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    # a tmp copy of the page with esr_ohm=0.2 removed fails (X), naming the keyword
    assert "esr_ohm=0.2, " in block
    stripped = tmp_path / "page.md"
    stripped.write_text(text.replace("esr_ohm=0.2, ", "", 1))
    stripped_reason = _block_defect(extract_tagged_block(stripped.read_text(), ESR_TAG))
    _p0(f"[WF-17 s14] (c-page) tags removed / doubled raise LookupError; esr_ohm removed -> "
        f"{stripped_reason}")
    assert stripped_reason is not None and "esr_ohm" in stripped_reason
    for token in (".json", ".npz", ".s4p", ".xdmf", "AED", "HFSS"):
        assert token not in text, token

    # (In): the block's inputs are the gate's; Study.solve patched (no mesh)
    calls = []

    def fake_solve(self, *args, **kwargs):
        calls.append(self)
        return SimpleNamespace(port_definitions=[], capacitor_impedance_ohm=None)

    monkeypatch.setattr(Study, "solve", fake_solve)
    t0 = time.perf_counter()
    ns, _ = exec_esr_block(monkeypatch)
    st = ns["study"]
    want_exc = CapacitorSheets(C_F, ESR_OHM, on=(PORTS[PLANT_INDEX],))
    _p0(f"[WF-17 s14] (In) solve calls {len(calls)}; geometry {st.geometry!r}; material "
        f"{st.material!r}; excitation {st.excitation!r}; frequency "
        f"{float(st.frequency_hz).hex()} (want {FREQUENCY_64_HZ.hex()}); degree {st.degree}; "
        f"conductor {st.conductor!r}; c_f hex {float(st.excitation.c_f).hex()} (C_F "
        f"{C_F.hex()}); exec {time.perf_counter() - t0:.2f} s")
    assert len(calls) == 1 and calls[0] is st
    assert type(st.geometry) is FourPortBirdcageGeometry and st.geometry == _geometry()
    assert st.geometry.strip_fraction == GATED_WIDTH_FRACTION
    assert st.material == _phantom()
    assert float(st.material.sigma).hex() == SALINE_SIGMA.hex()
    assert float(st.material.epsilon_r).hex() == SALINE_EPSILON_R.hex()
    assert type(st.excitation) is CapacitorSheets and st.excitation == want_exc
    assert float(st.excitation.c_f).hex() == C_F.hex()
    assert float(st.excitation.esr_ohm).hex() == ESR_OHM.hex()
    assert st.excitation.on == ("P2",) and PORTS[PLANT_INDEX] == "P2"
    assert st.excitation.drive == LumpedPortDrive()
    assert float(st.frequency_hz).hex() == FREQUENCY_64_HZ.hex()
    assert st.degree == 1
    assert type(st.conductor) is LeontovichWall and st.conductor == LeontovichWall(COPPER_SIGMA)
    assert float(st.conductor.sigma_s_per_m).hex() == COPPER_SIGMA.hex()


def test_wf17_s14_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s14] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _bands_defect(block) is None, _bands_defect(block)
    altered = block.replace(f"{REDUCTION_BAND:.0e}", f"{REDUCTION_BAND * 10:.0e}", 1)
    reason = _bands_defect(altered)
    _p0(f"[WF-17 s14] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "circuit reduction" in reason


def test_wf17_s14_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s14] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s14_esr_block_reproduces_wf16_s10(monkeypatch):
    comm = _solve_preamble()
    real_build = FourPortBirdcageGeometry.build
    cache = {"key": None, "built": None, "hits": 0, "real": 0}

    def _build(self, comm_):
        if cache["built"] is not None and cache["key"] == self and comm_ is comm:
            cache["hits"] += 1
            return cache["built"]
        cache["real"] += 1
        out = real_build(self, comm_)
        cache["key"], cache["built"] = self, out
        return out

    t_node = time.perf_counter()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        t0 = time.perf_counter()
        s0_res = _study(LumpedPortDrive()).solve(comm)
        t_s0 = time.perf_counter() - t0
        t0 = time.perf_counter()
        ns, _ = exec_esr_block(monkeypatch)
        t_page = time.perf_counter() - t0
        t0 = time.perf_counter()
        plant = _study(CapacitorSheets(C_F, PLANT_OHM, on=(PORTS[PLANT_INDEX],))).solve(comm)
        t_plant = time.perf_counter() - t0
    _p0(f"[WF-17 s14] build cache: {cache['real']} real build, {cache['hits']} hits (3 Studies); "
        f"S0 {t_s0:.1f} s, page block {t_page:.1f} s, S_plant {t_plant:.1f} s")

    # (R0)
    s0 = s0_res.s_matrix
    means = _class_means(s0)
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]
    rel_r0 = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    _p0(f"[WF-17 s14] (R0) cells {s0_res.cells} (record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); "
        "Cu class means vs record: " + " / ".join(f"{c} {rel_r0[c]:.3e}" for c in CLASSES)
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}; W1:976 {W1_R0})")

    # (I) the page's identity
    result = ns["result"]
    ids = [p.port_id for p in result.port_definitions]
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    pred = reduce_terminated_ports(s0, z0, {PLANT_INDEX: _rc(ESR_OHM)})
    r_i = _rel(result.s_matrix, pred)
    _p0(f"[WF-17 s14] (I) page ports {ids}; Z {result.capacitor_impedance_ohm!r} (want "
        f"{_rc(ESR_OHM)!r}); S {result.s_matrix.shape}; rel {r_i:.6e} "
        f"({r_i / REDUCTION_BAND:.4f} x REDUCTION_BAND {REDUCTION_BAND:.0e}; W1:980 {W1_I_P2})")

    # (N1) planted 0.4 ohm on P2
    r_n1 = _rel(plant.s_matrix, pred)
    _p0(f"[WF-17 s14] (N1) planted r = {PLANT_OHM:g} on P{PLANT_INDEX + 1}: {r_n1:.6e} = "
        f"{r_n1 / REDUCTION_BAND:.3f} x REDUCTION_BAND (ASSERTED > {REDUCTION_BAND:.0e}; "
        f"W1:996 {W1_N1})")
    _p0(f"[WF-17 s14] W-esr node {time.perf_counter() - t_node:.1f} s")

    assert cache["real"] == 1 and cache["hits"] == 2
    assert s0_res.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert max(rel_r0.values()) <= LEG_D0_REPRODUCTION_BAND
    assert ids == ["P1", "P3", "P4"]
    assert result.capacitor_impedance_ohm == _rc(ESR_OHM)
    assert r_i <= REDUCTION_BAND
    assert r_n1 > REDUCTION_BAND
