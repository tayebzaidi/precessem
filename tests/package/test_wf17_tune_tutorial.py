"""`WF-17` step 16: the tuning tutorial page reads ``Study.tune(k=1)``'s
``C_mode`` on the F-small PEC hole four-port from Python and from a ``tune:``
study file, and reproduces `WF-16` step 6(b)'s tuned capacitance, its lossless
Route-A circuit-reduction identity and its k = 2 control, and `WF-1` step 13's
study-file ``C_mode`` record.

Brief: docs/planning/chunks/WF-17.md, "§9 items 627, 632, 633, 634" (item 632),
over the template "§9 item 621".  The solve node runs only with
``FEM_EM_WF17_S16=1`` (and then asserts the complex build and ``-n
ANS6_EXACT_STRIP_RECORD_RANK_WIDTH``, never skips).  The other nodes build no
mesh.

Anchors: (X) one block per tag (three tags); the python block holds one
``Study(`` call, one ``.tune(`` call at k = 1, and no ``solve`` /
``run_study_file`` / ``CapacitorSheets`` call; the yaml block one
``run_study_file(`` call whose argument is the tracked file's repository path
and no ``Study`` / ``solve`` / ``tune`` call (AST names); (In) the python
block's Study and the file's parsed Study equal the gate's inputs exactly (no
mesh); (E) the bands block = the imported constants at ``.0e``; (R0) S0's
cells and class means vs the ``"pec"`` hole record at
``LEG_D0_REPRODUCTION_BAND``; (T) the python block's ``C_mode`` vs S0's own
``tune(k=1)`` and vs `TH-17`'s record (the gate's comparison) at
``REDUCTION_BAND``; (I-A) Route A at the block's ``C_mode`` (lossless, P2-P4)
vs the circuit reduction of S0 at ``REDUCTION_BAND``; (TY) the yaml block's
``summary["tune"]`` ``C_mode`` vs `WF-1` step 13's record
(``ITEM429_C_MODE_F``, imported by name from
``tests/package/test_wf1_yaml_capacitors.py``, compared as that module's (a')
does) at ``REDUCTION_BAND``; (D) 0 dead links.  Controls (rule (e)): (c1) the
k = 2 capacitance leaves ``|Im z_1|/|z_1|`` above ``REDUCTION_BAND``
(asserted, backed by the gate's 3.2883e-01); (c-page) by construction.  Bands
only are asserted; every reading is printed.  Route A and (c1) run here, not
on the page (the gap-3 ruling gives the lossless capacitor route no page).

Mesh economy (disclosed, as the gate's ``s6b`` fixture): inside the solve node
``FourPortBirdcageGeometry.build`` is wrapped so the first call builds the real
mesh and later equal-geometry calls on the same communicator return it.  The
yaml block names the file by its repository path; it is exec'd from the
scratch, so rank 0 copies the tracked file to the same relative path under the
scratch first (disclosed adaptation, as step 15; the copy's sha256 is printed
beside the tracked one).
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

import fem_em_solver
from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.ports.circuit import reduce_terminated_ports, s_to_z
from fem_em_solver.study import GATED_WIDTH_FRACTION
from fem_em_solver.study_yaml import _load, load_study

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_study_copper import _geometry
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
from tests.package.test_wf16_study_tune_capacitors import _phantom, _study
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf1_yaml_capacitors import ITEM429_C_MODE_F
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
from tests.validation.test_th17_birdcage_eigenmodes import _modal_z

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "tuning_a_study.md"
PAGE_ID = "docs/guide/tuning_a_study.md"
FILE_REL = "tests/package/wf1/f_small_hole_pec_tune.yaml"
YAML_FILE = REPO / FILE_REL
OTHER_REL = "tests/package/wf1/f_small_hole_pec_control.yaml"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s16"
S16_ENV = "FEM_EM_WF17_S16"
MIN_GUIDE_PAGES = 65
SHEETS = ("P2", "P3", "P4")

PY_TAG = "tune-python"
YAML_TAG = "tune-yaml"
BANDS_TAG = "tune-bands"

# Printed beside the readings (`20261015T201110Z_WF-16-s6b-W1.log.gz` :940 /
# :954 / :970); never asserted.
W1_C1 = "6.626796421e-11"
W1_C2 = "9.873468900e-11"
W1_IA = "1.391811e-05"
W1_C1_CONTROL = "3.2883e-01"


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


def _py_defect(block: str):
    """None when the python block has the (X) shape, else the reason."""
    n_study = len(_calls_named(block, "Study"))
    if n_study != 1:
        return f"{n_study} Study calls; need exactly one"
    tunes = [c for c in _calls_named(block, "tune") if isinstance(c.func, ast.Attribute)]
    if len(tunes) != 1:
        return f"{len(tunes)} .tune( calls; need exactly one"
    t = tunes[0]
    ks = [ast.literal_eval(a) for a in t.args] + [ast.literal_eval(k.value) for k in t.keywords
                                                  if k.arg == "k"]
    if ks != [1]:
        return f".tune( k arguments {ks!r} != [1]"
    for name in ("solve", "solve_sweep", "run_study_file", "CapacitorSheets"):
        if _calls_named(block, name):
            return f"the block calls {name}"
    return None


def _yaml_defect(block: str):
    """None when the study-file block has the (X) shape, else the reason."""
    runs = _calls_named(block, "run_study_file")
    if len(runs) != 1:
        return f"{len(runs)} run_study_file calls; need exactly one"
    arg = ast.literal_eval(runs[0].args[0])
    if arg != FILE_REL:
        return f"run_study_file argument {arg!r} != {FILE_REL!r}"
    for name in ("Study", "solve", "solve_sweep", "tune"):
        if _calls_named(block, name):
            return f"the block calls {name}"
    return None


def _bands_lines() -> list:
    return [
        f"Python C_mode vs the gate's tuned capacitance: {REDUCTION_BAND:.0e}",
        f"Study-file C_mode vs the study-file record: {REDUCTION_BAND:.0e}",
        f"Untuned 4×4 vs the perfect-conductor record: {LEG_D0_REPRODUCTION_BAND:.0e}",
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
    raw = os.environ.get(S16_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S16_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(f"{S16_ENV} unset: the tutorial solve window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port tune Study needs the complex build"
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def _prepare_scratch():
    """Rank 0 clears the scratch and copies the tracked file to its relative path."""
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    data = YAML_FILE.read_bytes()
    if comm.rank == 0:
        print(f"[WF-17 s16] page sha256 {_sha(raw)}", flush=True)
        if SCRATCH_ROOT.exists():
            shutil.rmtree(SCRATCH_ROOT)
        target = SCRATCH_ROOT / FILE_REL
        target.parent.mkdir(parents=True)
        target.write_bytes(data)
        print(f"[WF-17 s16] copied {FILE_REL} sha256 {_sha(target.read_bytes())} (tracked "
              f"{_sha(data)})", flush=True)
    comm.Barrier()


def exec_block(monkeypatch, tag):
    """chdir into the scratch on every rank, exec the tagged block as __main__."""
    block = extract_tagged_block(PAGE.read_text(), tag)
    _p0(f"[WF-17 s16] block {tag} sha256 {_sha(block.encode())}")
    monkeypatch.chdir(SCRATCH_ROOT)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{tag}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s16] block {tag} stdout (rank 0): {out!r}")
    return namespace, out


def _assert_gate_inputs(st, label):
    _p0(f"[WF-17 s16] (In) {label}: geometry {st.geometry!r}; material {st.material!r}; "
        f"excitation {st.excitation!r}; frequency {float(st.frequency_hz).hex()} (want "
        f"{FREQUENCY_64_HZ.hex()}); degree {st.degree}; conductor {st.conductor!r}")
    assert type(st.geometry) is FourPortBirdcageGeometry and st.geometry == _geometry()
    assert st.geometry.strip_fraction == GATED_WIDTH_FRACTION
    assert st.material == _phantom()
    assert float(st.material.sigma).hex() == SALINE_SIGMA.hex()
    assert float(st.material.epsilon_r).hex() == SALINE_EPSILON_R.hex()
    assert type(st.excitation) is LumpedPortDrive and st.excitation == LumpedPortDrive()
    assert float(st.frequency_hz).hex() == FREQUENCY_64_HZ.hex()
    assert st.degree == 1
    assert st.conductor is None


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s16_blocks_and_inputs(monkeypatch, tmp_path):
    text = _text()
    pblock = extract_tagged_block(text, PY_TAG)
    yblock = extract_tagged_block(text, YAML_TAG)
    extract_tagged_block(text, BANDS_TAG, lang="text")
    p_reason, y_reason = _py_defect(pblock), _yaml_defect(yblock)
    _p0(f"[WF-17 s16] (X) {PY_TAG}: Study {len(_calls_named(pblock, 'Study'))}, tune "
        f"{len(_calls_named(pblock, 'tune'))}, solve {len(_calls_named(pblock, 'solve'))}, "
        f"run_study_file {len(_calls_named(pblock, 'run_study_file'))}; defect {p_reason}")
    _p0(f"[WF-17 s16] (X) {YAML_TAG}: run_study_file "
        f"{len(_calls_named(yblock, 'run_study_file'))}, Study {len(_calls_named(yblock, 'Study'))}, "
        f"tune {len(_calls_named(yblock, 'tune'))}; defect {y_reason}")
    assert p_reason is None, p_reason
    assert y_reason is None, y_reason
    assert pblock.count("Study(") == 1 and "run_study_file" not in pblock
    assert "Study(" not in yblock and yblock.count("run_study_file(") == 1
    assert "(capacitor_sheets_with_esr.md)" in text
    assert "(tuning_from_yaml.md)" in text
    assert f"(../../{FILE_REL})" in text
    assert "lossless" in text

    # (c-page), by construction: a removed / doubled tag raises
    for tag, lang in ((PY_TAG, "python"), (YAML_TAG, "python"), (BANDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    # a tmp copy of the page tuning k = 2 fails (X); one naming another file fails (X)
    assert "study.tune(1)" in pblock
    k2 = tmp_path / "page.md"
    k2.write_text(text.replace("study.tune(1)", "study.tune(2)", 1))
    k2_reason = _py_defect(extract_tagged_block(k2.read_text(), PY_TAG))
    other = tmp_path / "page2.md"
    other.write_text(text.replace(f'run_study_file("{FILE_REL}")',
                                  f'run_study_file("{OTHER_REL}")', 1))
    other_reason = _yaml_defect(extract_tagged_block(other.read_text(), YAML_TAG))
    _p0(f"[WF-17 s16] (c-page) tags removed / doubled raise LookupError for all three blocks; "
        f"k = 2 -> {k2_reason}; other file -> {other_reason}")
    assert k2_reason is not None and "k arguments" in k2_reason
    assert other_reason is not None and "f_small_hole_pec_control" in other_reason
    for token in (".json", ".npz", ".s4p", ".xdmf", "AED", "HFSS"):
        assert token not in text, token

    # (In) the python block: Study.tune patched (no mesh), Study.solve must not run
    calls = []

    def fake_tune(self, *args, **kwargs):
        calls.append((self, args, kwargs))
        return SimpleNamespace(k=1, c_mode_farad=0.0)

    def no_solve(self, *args, **kwargs):
        raise AssertionError("Study.solve must not run in (In)")

    monkeypatch.setattr(Study, "tune", fake_tune)
    monkeypatch.setattr(Study, "solve", no_solve)
    _prepare_scratch()
    t0 = time.perf_counter()
    ns, _ = exec_block(monkeypatch, PY_TAG)
    st = ns["study"]
    _p0(f"[WF-17 s16] (In) tune calls {[(a, k) for _, a, k in calls]}; exec "
        f"{time.perf_counter() - t0:.2f} s")
    assert len(calls) == 1 and calls[0][0] is st and calls[0][1] == (1,) and calls[0][2] == {}
    _assert_gate_inputs(st, PY_TAG)

    # (In) the study file: the same Study, tune k = 1 (parse only, no mesh)
    sy = load_study(YAML_FILE)
    _assert_gate_inputs(sy, FILE_REL)
    outputs = _load(YAML_FILE)[2]
    _p0(f"[WF-17 s16] (In) {FILE_REL} outputs {outputs!r}")
    assert outputs.get("tune") == {"k": 1}


def test_wf17_s16_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s16] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _bands_defect(block) is None, _bands_defect(block)
    altered = block.replace(f"{LEG_D0_REPRODUCTION_BAND:.0e}",
                            f"{LEG_D0_REPRODUCTION_BAND * 10:.0e}", 1)
    reason = _bands_defect(altered)
    _p0(f"[WF-17 s16] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "perfect-conductor record" in reason


def test_wf17_s16_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s16] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s16_tune_blocks_reproduce_wf16_s6b(monkeypatch):
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

    calls = {"n": 0, "paths": []}
    original = fem_em_solver.run_study_file

    def spy(path, *args, **kwargs):
        calls["n"] += 1
        calls["paths"].append(str(path))
        return original(path, *args, **kwargs)

    t_node = time.perf_counter()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        mp.setattr(fem_em_solver, "run_study_file", spy)
        t0 = time.perf_counter()
        s0 = _study(LumpedPortDrive())
        s0_res = s0.solve(comm)
        t1 = s0.tune(k=1, comm=comm)
        t2 = s0.tune(k=2, comm=comm)
        assert s0.result is s0_res
        t_s0 = time.perf_counter() - t0
        _prepare_scratch()
        t0 = time.perf_counter()
        ns_p, _ = exec_block(monkeypatch, PY_TAG)
        t_py = time.perf_counter() - t0
        c_page = float(ns_p["tuned"].c_mode_farad)
        exc_a = CapacitorSheets(c_page, 0.0, on=SHEETS)
        t0 = time.perf_counter()
        res_a = _study(exc_a).solve(comm)
        t_a = time.perf_counter() - t0
        t0 = time.perf_counter()
        ns_y, _ = exec_block(monkeypatch, YAML_TAG)
        t_yaml = time.perf_counter() - t0
    _p0(f"[WF-17 s16] build cache: {cache['real']} real build, {cache['hits']} hits (4 Studies); "
        f"S0 + two tunes {t_s0:.1f} s, python block {t_py:.1f} s, Route A {t_a:.1f} s, "
        f"yaml block {t_yaml:.1f} s")

    # (R0)
    s0m = np.asarray(s0_res.s_matrix, dtype=np.complex128)
    means = _class_means(s0m)
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]
    rel_r0 = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    _p0(f"[WF-17 s16] (R0) cells {s0_res.cells} (record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); "
        "PEC class means vs record: " + " / ".join(f"{c} {rel_r0[c]:.3e}" for c in CLASSES)
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})")

    # (T) the python block's C_mode
    tuned = ns_p["tuned"]
    c1 = float(t1.c_mode_farad)
    r_t = abs(c_page / c1 - 1.0)
    a_t = c_page / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(f"[WF-17 s16] (T) block k {tuned.k}, C_mode {c_page:.9e} F (hex {c_page.hex()}); S0 "
        f"tune C_1 {c1:.9e} F (hex {c1.hex()}): rel {r_t:.3e} (ASSERTED <= {REDUCTION_BAND:.0e}); "
        f"vs TH-17 record {a_t:+.6e} (ASSERTED |.| <= {REDUCTION_BAND:.0e}); gate W1:940 C_1 "
        f"{W1_C1}, C_2 {W1_C2}; S0 C_2 {float(t2.c_mode_farad):.9e} F (printed)")

    # (I-A) Route A at the block's C_mode, lossless, P2-P4
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    zc = exc_a.impedance_ohm(FREQUENCY_64_HZ)
    pa = reduce_terminated_ports(s0m, z0, {1: zc, 2: zc, 3: zc})
    ids_a = [p.port_id for p in res_a.port_definitions]
    r_ia = abs(res_a.s_matrix[0, 0] - pa[0, 0]) / abs(pa[0, 0])
    _p0(f"[WF-17 s16] (I-A) ports {ids_a}; Z {res_a.capacitor_impedance_ohm!r} (want {zc!r}); "
        f"S11 {complex(res_a.s_matrix[0, 0]):.9e} vs reduction {complex(pa[0, 0]):.9e}: rel "
        f"{r_ia:.6e} ({r_ia / REDUCTION_BAND:.4f} x REDUCTION_BAND {REDUCTION_BAND:.0e}; "
        f"W1:954 {W1_IA})")

    # (c1) the k = 2 capacitance leaves Im z_1 of order 0.33 |z_1|
    zc2 = CapacitorSheets(float(t2.c_mode_farad), 0.0, on=SHEETS).impedance_ohm(FREQUENCY_64_HZ)
    z1 = complex(t1.z_k[1])
    zt2 = _modal_z(s_to_z(s0m, z0) + zc2 * np.eye(4))[0]
    r_c1 = abs(zt2[1].imag) / abs(z1)
    _p0(f"[WF-17 s16] (c1) |Im z_1(C_2)|/|z_1| = {r_c1:.4e} ({r_c1 / REDUCTION_BAND:.1f} x "
        f"REDUCTION_BAND; ASSERTED > {REDUCTION_BAND:.0e}; W1:970 {W1_C1_CONTROL})")

    # (TY) the yaml block's summary["tune"]
    run = ns_y["run"]
    tsum = run.summary["tune"]
    c_y = float(tsum["c_mode_farad"])
    a_y = ITEM429_C_MODE_F / c_y - 1.0
    _p0(f"[WF-17 s16] (TY) spy calls {calls['n']} {calls['paths']}; family "
        f"{run.summary['family']!r}; cells {run.summary['cells']}; tune k {tsum['k']}, C_mode "
        f"{c_y:.9e} F (hex {c_y.hex()}); ITEM429_C_MODE_F {ITEM429_C_MODE_F.hex()} vs in-run: "
        f"{a_y:+.3e} (ASSERTED |.| <= {REDUCTION_BAND:.0e}); vs python block "
        f"{c_y / c_page - 1.0:+.3e}, hex equal {c_y.hex() == c_page.hex()} (printed)")
    _p0(f"[WF-17 s16] W-tune node {time.perf_counter() - t_node:.1f} s")

    assert cache["real"] == 1 and cache["hits"] == 3
    assert calls["n"] == 1 and calls["paths"] == [FILE_REL]
    assert s0_res.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert max(rel_r0.values()) <= LEG_D0_REPRODUCTION_BAND
    assert tuned.k == 1 and t1.k == 1 and t2.k == 2
    assert r_t <= REDUCTION_BAND
    assert abs(a_t) <= REDUCTION_BAND
    assert ids_a == ["P1"] and res_a.s_matrix.shape == (1, 1)
    assert res_a.capacitor_impedance_ohm == zc
    assert r_ia <= REDUCTION_BAND
    assert zc2 != zc
    assert r_c1 > REDUCTION_BAND
    assert run.summary["family"] == "four_port_birdcage"
    assert run.summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert tsum["k"] == 1
    assert abs(a_y) <= REDUCTION_BAND
