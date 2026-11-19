"""`WF-17` step 7: the thermal-box tutorial page runs `WF-1`'s mri:29 /
mri:31 / mri:32 file sets through ``run_study_file`` and reproduces
`MAT-21`'s P36 / P38 records.

Brief: docs/planning/chunks/WF-17.md, "§9 item 563" and the review ruling
2026-10-20 03:00.  The three solve nodes run only with ``FEM_EM_WF17_S7=1``
(and then assert the complex build and ``-n 2``, never skip).  The other nodes
build no mesh.

Anchors: (X) one block per tag, one ``run_study_file(`` call, no ``Study(``;
(F) each block's ``FILES`` = its example's docstring file set; (M) the
malformed files refuse naming their field; (E) the records block = P36 / P38's
digits; (R-run) each block runs each file once; (P36) / (P38) at
``INTERIOR_BAND`` and the records' ``.6e`` digits, Gauss at ``S2_GAUSS_BAND``,
balance at ``E_BAND``; (D) 0 dead links.  Controls (rule (e)): (c-src) and
(c-step) asserted, backed; (c-page) by construction.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from mpi4py import MPI

import fem_em_solver
from fem_em_solver import cli, load_study
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_yaml_thermal_ball_probes import MALFORMED as MALFORMED_BALL_PROBE
from tests.package.test_wf1_yaml_thermal_box import MALFORMED_DUTY
from tests.package.test_wf1_yaml_thermal_box_probes import MALFORMED as MALFORMED_BOX_PROBE
from tests.package.test_wf1_yaml_thermal_box_quadratic import P36_C_SRC, P36_S1, P36_T, _field
from tests.package.test_wf1_yaml_thermal_box_quadratic_duty import P38_C_STEP, P38_T_QDUTY
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND
from tests.validation.test_mat21_box_fem import T_STAR, _interior_chi, _steady_ex
from tests.validation.test_mat21_box_fem_p2 import _exact_p, _rel_p
from tests.validation.test_mat21_box_quadratic_duty_closed_form import _quad_duty
from tests.validation.test_mat21_box_quadratic_duty_fem import SCHEDULE
from tests.validation.test_mat21_box_quadratic_fem import _quad_ex, _quad_trans_ex

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "thermal_box_from_yaml.md"
PAGE_ID = "docs/guide/thermal_box_from_yaml.md"
WF1_DIR = REPO / "tests" / "package" / "wf1"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s7"
S7_ENV = "FEM_EM_WF17_S7"
MIN_GUIDE_PAGES = 52

EXAMPLES = {
    "thermal-box-mri29": REPO / "examples" / "mri" / "29_femem_run_f2182_gel_box_files.py",
    "thermal-box-mri31": REPO
    / "examples"
    / "mri"
    / "31_femem_run_gel_box_source_and_probe_files.py",
    "thermal-box-mri32": REPO
    / "examples"
    / "mri"
    / "32_femem_run_ball_probes_and_duty_quadratic_files.py",
}
CHECK_TAG = "thermal-box-check"
RECORDS_TAG = "thermal-box-records"
QUAD_STEADY = "thermal_box_quadratic_steady_m4.yaml"
QUAD_CN = "thermal_box_quadratic_cn_m4.yaml"
QUAD_DUTY = "thermal_box_quadratic_duty_m4.yaml"
N_SOLVE_FILES = 12


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _files_of(block: str) -> list:
    """``FILES`` from the parsed AST (never by exec)."""
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "FILES" for t in node.targets
        ):
            return list(ast.literal_eval(node.value))
    raise LookupError("no FILES assignment in the block")


def _docstring_files(example: Path) -> list:
    doc = ast.get_docstring(ast.parse(example.read_text()))
    return list(dict.fromkeys(re.findall(r"[\w]+\.yaml", doc)))


def _solve_files(example: Path) -> list:
    return [f for f in _docstring_files(example) if not f.startswith("malformed_")]


def _calls_named(block: str, name: str) -> int:
    n = 0
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Call):
            fn = node.func
            fname = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)
            if fname == name:
                n += 1
    return n


def _file_set_defect(files: list, example: Path):
    """None when ``files`` is the example's solve-file set, else the reason."""
    want = _solve_files(example)
    if len(files) != len(set(files)):
        return f"duplicates in FILES {files}"
    missing = [f for f in want if f not in files]
    extra = [f for f in files if f not in want]
    if missing or extra:
        return f"missing {missing}, extra {extra}"
    return None


def _knob() -> bool:
    import os

    raw = os.environ.get(S7_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S7_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(
            f"{S7_ENV} unset: the tutorial solve windows are environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "P36 / P38 are complex-build records"
    assert comm.size == 2, comm.size
    return comm


def exec_thermal_block(tag: str, monkeypatch):
    """Copy the block's files into the scratch, chdir there, exec the block.

    Returns ``(namespace, calls, stdout)``; the spy calls through.
    """
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), tag)
    files = _files_of(block)
    scratch = SCRATCH_ROOT / tag
    if comm.rank == 0:
        print(
            f"[WF-17 s7] page sha256 {_sha(raw)}; block {tag} sha256 {_sha(block.encode())}",
            flush=True,
        )
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        for name in files:
            data = (WF1_DIR / name).read_bytes()
            (scratch / name).write_bytes(data)
            print(f"[WF-17 s7] copied {name} sha256 {_sha(data)}", flush=True)
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
    _p0(f"[WF-17 s7] block {tag} stdout (rank 0): {out!r}")
    return namespace, calls, out, files


def _check_r_run(namespace, calls, files):
    runs = namespace["runs"]
    fam = {name: runs[name].summary["family"] for name in files}
    _p0(
        f"[WF-17 s7] (R-run) runs {len(runs)} / spy calls {calls['n']} / FILES {len(files)}; "
        f"families {fam}"
    )
    assert len(runs) == len(files)
    assert calls["n"] == len(files)
    for name in files:
        want = "thermal_ball" if name.startswith("thermal_ball") else "thermal_box"
        assert fam[name] == want, (name, fam[name])
    return runs


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s7_blocks_and_file_sets():
    text = _text()
    union = []
    for tag, example in EXAMPLES.items():
        block = extract_tagged_block(text, tag)
        n_run = _calls_named(block, "run_study_file")
        n_study = _calls_named(block, "Study")
        files = _files_of(block)
        defect = _file_set_defect(files, example)
        _p0(
            f"[WF-17 s7] (X/F) {tag}: run_study_file calls {n_run}, Study calls {n_study}; "
            f"FILES {files}; docstring set {_solve_files(example)}; defect {defect}"
        )
        assert block.count("run_study_file(") == 1
        assert n_run == 1 and n_study == 0 and "Study(" not in block
        assert defect is None, defect
        for name in files:
            assert (WF1_DIR / name).is_file(), name
        union += files
    assert len(union) == N_SOLVE_FILES == len(set(union))
    extract_tagged_block(text, CHECK_TAG, lang="bash")
    extract_tagged_block(text, RECORDS_TAG, lang="text")
    # (c-page), by construction: a removed / doubled tag, a dropped file
    tag0 = "thermal-box-mri31"
    with pytest.raises(LookupError):
        extract_tagged_block(text.replace(f"```python {tag0}", "```python other"), tag0)
    opener = f"```python {tag0}\n"
    doubled = text + "\n" + opener + extract_tagged_block(text, tag0) + "```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, tag0)
    dropped = _files_of(extract_tagged_block(text, tag0))[1:]
    reason = _file_set_defect(dropped, EXAMPLES[tag0])
    _p0(f"[WF-17 s7] (c-page) one file dropped -> {reason}")
    assert reason is not None and QUAD_STEADY in reason
    for name in ("Study(", ".json", ".xdmf", ".npz"):
        assert name not in text.replace("run_study_file(", ""), name


def test_wf17_s7_malformed_files_refuse(capsys):
    comm = MPI.COMM_WORLD
    check = extract_tagged_block(_text(), CHECK_TAG, lang="bash")
    cases = [
        (MALFORMED_DUTY, "thermal.duty"),
        (MALFORMED_BOX_PROBE, "outputs.probes[1].radius"),
        (MALFORMED_BALL_PROBE, "outputs.probes[1].radius"),
    ]
    for path, field in cases:
        assert f"femem check {path.name}" in check, path.name
        got = _field(lambda: load_study(path))
        rc = cli.main(["check", str(path)])
        rcs = comm.allgather(rc)
        _p0(f"[WF-17 s7] (M) {path.name}: field {got!r} (want {field!r}); rc per rank {rcs}")
        assert got == field
        assert all(r == 2 for r in rcs)


def test_wf17_s7_records_on_the_page():
    block = extract_tagged_block(_text(), RECORDS_TAG, lang="text")
    values = [line.rsplit("= ", 1)[1].strip() for line in block.splitlines() if "= " in line]
    _p0(f"[WF-17 s7] (E) records block {values} (P36 {P36_S1} / {P36_T}; P38 {P38_T_QDUTY})")
    assert values == [P36_S1, P36_T, P38_T_QDUTY]


def test_wf17_s7_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s7] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve windows (-n 2)


def test_wf17_s7_mri29_block_runs(monkeypatch):
    _solve_preamble()
    t0 = time.perf_counter()
    ns, calls, _, files = exec_thermal_block("thermal-box-mri29", monkeypatch)
    runs = _check_r_run(ns, calls, files)
    for name in files:
        st = runs[name].summary["temperature"]
        keys = {k: st[k] for k in ("source", "duty", "max_balance_residual") if k in st}
        _p0(f"[WF-17 s7] mri:29 {name} summary (printed beside, MAT-21 step 1 / P26's): {keys}")
    _p0(f"[WF-17 s7] W29 leg {time.perf_counter() - t0:.1f}s")


def test_wf17_s7_mri31_block_reproduces_p36(monkeypatch):
    comm = _solve_preamble()
    t0 = time.perf_counter()
    ns, calls, _, files = exec_thermal_block("thermal-box-mri31", monkeypatch)
    runs = _check_r_run(ns, calls, files)
    for name in files:
        if "probes" in name:
            _p0(
                f"[WF-17 s7] mri:31 {name} probes (real-build records, printed beside only): "
                f"{runs[name].summary.get('probes')}"
            )
    # steady (the gate's read: the submesh)
    run = runs[QUAD_STEADY]
    th = run.thermal
    msh = th.submesh
    th_q = th.temperature - T_INF
    chi = _interior_chi(msh)
    ex_q, npts = _exact_p(msh, _quad_ex())
    s1 = _rel_p(msh, th_q, ex_q, chi)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=th.facet_tags)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_q * ds)))
    sp = float(np.real(th.source_power))
    gauss = abs(conv - sp) / sp
    ex_u, _ = _exact_p(msh, _steady_ex())
    c_src = _rel_p(msh, th_q, ex_u, chi)
    # CN (the gate's read: the mesh)
    run_t = runs[QUAD_CN]
    tht = run_t.thermal
    msh_t = tht.mesh
    ex_t, _ = _exact_p(msh_t, _quad_trans_ex(T_STAR))
    t_int = _rel_p(msh_t, tht.temperature - T_INF, ex_t, _interior_chi(msh_t))
    bal = run_t.summary["temperature"]["max_balance_residual"]
    _p0(
        f"[WF-17 s7] (P36 S1-quad) steady interior {s1:.6e} (P36 {P36_S1}; band "
        f"{INTERIOR_BAND:.0e}); {npts} P3 points; source "
        f"{run.summary['temperature']['source']!r}"
    )
    _p0(f"[WF-17 s7] (P36 S2-quad) Gauss {gauss:.3e} (band {S2_GAUSS_BAND:.0e})")
    _p0(
        f"[WF-17 s7] (P36 T-quad) CN interior at t* {t_int:.6e} (P36 {P36_T}; band "
        f"{INTERIOR_BAND:.0e})"
    )
    _p0(f"[WF-17 s7] (P36 E-quad) CN max balance residual {bal:.3e} (band {E_BAND:.0e})")
    _p0(
        f"[WF-17 s7] control (c-src) quadratic run vs uniform series {c_src:.4e} (P36 "
        f"{P36_C_SRC:.4e}; asserted > {INTERIOR_BAND:.0e})"
    )
    _p0(f"[WF-17 s7] W31 leg {time.perf_counter() - t0:.1f}s")
    assert run.summary["temperature"]["source"] == "transverse_quadratic"
    assert s1 <= INTERIOR_BAND
    assert f"{s1:.6e}" == P36_S1
    assert gauss <= S2_GAUSS_BAND
    assert t_int <= INTERIOR_BAND
    assert f"{t_int:.6e}" == P36_T
    assert bal <= E_BAND
    assert c_src > INTERIOR_BAND


def test_wf17_s7_mri32_block_reproduces_p38(monkeypatch):
    _solve_preamble()
    t0 = time.perf_counter()
    ns, calls, _, files = exec_thermal_block("thermal-box-mri32", monkeypatch)
    runs = _check_r_run(ns, calls, files)
    for name in files:
        if "probes" in name:
            _p0(
                f"[WF-17 s7] mri:32 {name} probes (real-build records, printed beside only): "
                f"{runs[name].summary.get('probes')}"
            )
    run = runs[QUAD_DUTY]
    th = run.thermal
    msh = th.mesh
    st = run.summary["temperature"]
    s0 = "".join(str(int(m)) for m in st["source_multipliers"])
    th_t = th.temperature - T_INF
    chi = _interior_chi(msh)
    ex_d, npts = _exact_p(msh, lambda x, y, z: _quad_duty(x, y, z, T_STAR))
    t_qd = _rel_p(msh, th_t, ex_d, chi)
    ex_s, _ = _exact_p(msh, _quad_trans_ex(T_STAR))
    c_step = _rel_p(msh, th_t, ex_s, chi)
    bal = st["max_balance_residual"]
    _p0(f"[WF-17 s7] (P38 S0) schedule {s0} (SCHEDULE {SCHEDULE})")
    _p0(
        f"[WF-17 s7] (P38 T-qduty) interior at t* {t_qd:.6e} (P38 {P38_T_QDUTY}; band "
        f"{INTERIOR_BAND:.0e}); {npts} P3 points"
    )
    _p0(f"[WF-17 s7] (P38 E-qduty) CN max balance residual {bal:.3e} (band {E_BAND:.0e})")
    _p0(
        f"[WF-17 s7] control (c-step) duty run vs unscheduled series {c_step:.4e} (P38 "
        f"{P38_C_STEP:.4e}; asserted > {INTERIOR_BAND:.0e})"
    )
    _p0(f"[WF-17 s7] W32 leg {time.perf_counter() - t0:.1f}s")
    assert s0 == SCHEDULE
    assert t_qd <= INTERIOR_BAND
    assert f"{t_qd:.6e}" == P38_T_QDUTY
    assert bal <= E_BAND
    assert c_step > INTERIOR_BAND
