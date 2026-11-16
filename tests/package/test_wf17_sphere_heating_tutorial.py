"""`WF-17` step 13: the gel-sphere heating page runs `WF-1`'s thermal and
probe files through ``run_study_file`` and reproduces `MAT-10`'s SAR-to-heat
record and `MAT-19`'s whole-ball probe identity.

Brief: docs/planning/chunks/WF-17.md, "§9 item 620".  The solve nodes run only
with ``FEM_EM_WF17_S13=heat`` (the heating node) or ``=probes`` (the probes
node); selected, a node asserts the complex build and ``-n 2`` and never skips.
The other nodes build no mesh.

Anchors: (X) one block per tag (five tags), one ``run_study_file(`` call and no
``Study(`` by AST name, the file named in each block, the bash line's tokens;
(E) the keys block = the three summary keys, the bands block = the four lines
rendered from the imported bands; (H) on the heating block's run, computed as
the gate ``test_wf1_yaml_thermal.py`` computes them: sub cells, (A) / (P) at
``.4e`` against the gate's own strings, (A) / (P) / (B) / (C) within their
bands; (Pr) on the probes block's run, as ``test_wf1_yaml_probes.py``: cells,
whole-ball value and volume within ``B_BAND``, summary values ==
``ball_average`` exactly; (D) 0 dead links.  Controls (rule (e)): (c-2G) the
control file's 2G cooling on the already-solved Study misses ``A_BAND``,
asserted, backed; (c-c2) the c2 ball 0.1 mm out misses by > ``WIDTH_BAND``,
asserted, backed; (c-page) by construction.  Every band, record, constant and
comparator is imported, none restated; readings are printed.
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
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

import fem_em_solver
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal.probe import ball_average
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf15_study_probe import C2_OUT
from tests.package.test_wf15_study_temperature import (
    MAT10_A_PRINT,
    MAT10_P_PRINT,
    MAT10_SUB_CELLS,
)
from tests.package.test_wf17_first_study_tutorial import LINK_PASS_RE, _checker_stdout, extract_tagged_block
from tests.package.test_wf1_yaml_probes import PROBES_FILE
from tests.package.test_wf1_yaml_thermal import CONTROL_2G_FILE, THERMAL_FILE, _cooling
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    A_BAND,
    B_BAND,
    C_BAND,
    G,
    H,
    K,
    P_BAND,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _theta_cf,
)
from tests.validation.test_mat19_virtual_temperature_probe import N_SUB, WIDTH_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "gel_sphere_heating_from_yaml.md"
PAGE_ID = "docs/guide/gel_sphere_heating_from_yaml.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s13"
S13_ENV = "FEM_EM_WF17_S13"
S13_VALUES = ("heat", "probes")
MIN_GUIDE_PAGES = 62

HEAT_TAG = "sphere-heating"
PROBES_TAG = "sphere-probes"
CLI_TAG = "sphere-heating-cli"
KEYS_TAG = "sphere-heating-keys"
BANDS_TAG = "sphere-heating-bands"
CLI_TOKENS = ["femem", "run", "sphere_mat10_thermal.yaml", "--output-dir",
              "femem_out/sphere_mat10_thermal"]
KEYS = ["source_power_w", "convected_power_w", "theta_max_k"]
BLOCK_FILES = {HEAT_TAG: THERMAL_FILE, PROBES_TAG: PROBES_FILE}

# Printed beside the readings, never asserted: `20261014T005614Z_WF-1-s10a-RG3-thermal.log`
# :149 / :150 and `20261020T201843Z_WF-1-s21-RGP-sphere-probes.log` :153 / :235.
RG3_R = "(A) 3.2463961454e-02, (P) 3.1923343273e-03, (B) 1.003e-15, (C) 1.107e-13"
RG3_C1 = "1.0649e+00"
RGP_A = "1.359e-15 / 1.221e-15"
RGP_C1 = "7.5333e-03"


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
        f"Interior error vs the closed form: {A_BAND:.1e}",
        f"Interior error, direct P1 comparator: {P_BAND:.0e}",
        f"Source power vs the EM dissipation: {B_BAND:.0e}",
        f"Whole-ball probe vs the volume mean: {B_BAND:.0e}",
    ]


def _lines_defect(block: str, want: list):
    """None when the block's lines equal ``want``, else the reason (naming the line)."""
    got = block.strip("\n").splitlines()
    if len(got) != len(want):
        missing = [w for w in want if w not in got]
        return f"{len(got)} lines != {len(want)}; missing {missing}"
    for g, w in zip(got, want):
        if g != w:
            return f"line {g!r} != {w!r}"
    return None


def _knob():
    raw = os.environ.get(S13_ENV, "").strip()
    if raw == "":
        return None
    if raw not in S13_VALUES:
        raise ValueError(f"{S13_ENV}={raw!r}: only unset/empty, 'heat' or 'probes' are accepted")
    return raw


def _solve_preamble(which: str):
    knob = _knob()
    if knob != which:
        pytest.skip(f"{S13_ENV}={knob!r}: the {which!r} tutorial solve window is "
                    "environment-selected (complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the gel-sphere EM solve needs the complex build"
    assert comm.size == 2, comm.size
    return comm


def exec_sphere_block(tag: str, monkeypatch):
    """Copy the block's file into the scratch, chdir there, exec the block.

    Returns ``(namespace, calls, stdout)``; the spy calls through.
    """
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), tag)
    name = _file_arg(block)
    src = BLOCK_FILES[tag]
    scratch = SCRATCH_ROOT / tag
    if comm.rank == 0:
        print(f"[WF-17 s13] page sha256 {_sha(raw)}; block {tag} sha256 {_sha(block.encode())}",
              flush=True)
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        data = (src.parent / name).read_bytes()
        (scratch / name).write_bytes(data)
        print(f"[WF-17 s13] copied {name} sha256 {_sha(data)} (tracked "
              f"{_sha(src.read_bytes())})", flush=True)
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
    _p0(f"[WF-17 s13] block {tag} stdout (rank 0): {out!r}")
    return namespace, calls, out


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s13_blocks_and_files():
    text = _text()
    for tag, path in BLOCK_FILES.items():
        block = extract_tagged_block(text, tag)
        n_run = len(_calls_named(block, "run_study_file"))
        n_study = len(_calls_named(block, "Study"))
        name = _file_arg(block)
        _p0(f"[WF-17 s13] (X) {tag}: run_study_file calls {n_run}, Study calls {n_study}; "
            f"file {name!r} (want {path.name!r})")
        assert block.count("run_study_file(") == 1
        assert n_run == 1 and n_study == 0 and "Study(" not in block
        assert name == path.name
        assert path.is_file()
        assert f"(../../tests/package/wf1/{path.name})" in text
    assert CONTROL_2G_FILE.is_file()
    cli = extract_tagged_block(text, CLI_TAG, lang="bash")
    cli_lines = [ln for ln in cli.splitlines() if ln.strip()]
    _p0(f"[WF-17 s13] (X) {CLI_TAG}: {cli_lines!r}")
    assert len(cli_lines) == 1
    assert cli_lines[0].split() == CLI_TOKENS
    extract_tagged_block(text, KEYS_TAG, lang="text")
    extract_tagged_block(text, BANDS_TAG, lang="text")
    # (c-page), by construction: a removed / doubled tag raises
    tags = ((HEAT_TAG, "python"), (PROBES_TAG, "python"), (CLI_TAG, "bash"),
            (KEYS_TAG, "text"), (BANDS_TAG, "text"))
    for tag, lang in tags:
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    _p0(f"[WF-17 s13] (c-page) removed / doubled tags raise LookupError for all {len(tags)} blocks")
    for token in ("Study(", ".json", ".npz", ".s4p", ".xdmf", ".h5", ".bp"):
        assert token not in text.replace("run_study_file(", ""), token


def test_wf17_s13_keys_and_bands_on_the_page():
    text = _text()
    keys = extract_tagged_block(text, KEYS_TAG, lang="text")
    bands = extract_tagged_block(text, BANDS_TAG, lang="text")
    _p0(f"[WF-17 s13] (E) keys block {keys.strip().splitlines()} (want {KEYS})")
    _p0(f"[WF-17 s13] (E) bands block {bands.strip().splitlines()} (rendered {_bands_lines()})")
    assert _lines_defect(keys, KEYS) is None, _lines_defect(keys, KEYS)
    assert _lines_defect(bands, _bands_lines()) is None, _lines_defect(bands, _bands_lines())
    # (c-page), by construction: one key dropped fails (E) and is named
    dropped = keys.replace("convected_power_w\n", "", 1)
    reason = _lines_defect(dropped, KEYS)
    _p0(f"[WF-17 s13] (c-page) one key dropped -> {reason}")
    assert reason is not None and "convected_power_w" in reason
    altered = bands.replace(f"{P_BAND:.0e}", f"{P_BAND * 10:.0e}", 1)
    reason_b = _lines_defect(altered, _bands_lines())
    _p0(f"[WF-17 s13] (c-page) one band exponent altered -> {reason_b}")
    assert reason_b is not None and "direct P1 comparator" in reason_b


def test_wf17_s13_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s13] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve windows (-n 2)


def test_wf17_s13_heating_block_reproduces_mat10(monkeypatch):
    comm = _solve_preamble("heat")
    t0 = time.perf_counter()
    ns, calls, _ = exec_sphere_block(HEAT_TAG, monkeypatch)
    t_block = time.perf_counter() - t0
    run = ns["run"]
    s = run.summary
    study, res, th = run.study, run.result, run.thermal
    n_calls, paths = calls["n"], list(calls["paths"])
    assert th is not None
    assert study.result is res, "(Once) a second EM solve happened"

    # (H) MAT-10's gates on the block's θ, computed as test_wf1_yaml_thermal.py does.
    T = th.temperature
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub, ft = th.submesh, th.facet_tags
    n_sub = sub.topology.index_map(3).size_global
    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    err_A = _interior_rel(sub, chi, (T - T_INF) - theta, theta)
    res_P = solve_steady_conduction(sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_P = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta, theta)
    int_q = float(np.real(_allsum(comm, th.heat_source * ufl.dx(domain=sub))))
    p_parent = th.cooling.source_scale * study.mean_sar()["dissipated_power_w"]
    err_B = abs(int_q - p_parent) / abs(p_parent)
    err_C = abs(th.convected_power - th.source_power) / abs(th.source_power)

    # (c-2G) the control file's cooling on the already-solved Study, no re-solve.
    c2 = _cooling(CONTROL_2G_FILE)
    th2 = study.temperature(c2)
    once_c = run.study.result is run.result
    sub2 = th2.submesh
    theta2 = _theta_cf(sub2, q_cf)
    c_2g = _interior_rel(sub2, _interior_indicator(sub2), (th2.temperature - T_INF) - theta2,
                         theta2)

    _p0(f"[WF-17 s13] (H-run) block spy calls {n_calls} {paths}; family {s['family']!r}; "
        f"degree {s['degree']!r}; temperature keys {list(s['temperature'])}; "
        f"summary['temperature'] {s['temperature']}")
    _p0(f"[WF-17 s13] (H) sub cells {n_sub} (record {MAT10_SUB_CELLS}); (A) = {err_A:.10e} "
        f"[{err_A:.4e}] vs {MAT10_A_PRINT} (band {A_BAND}); (P) = {err_P:.10e} [{err_P:.4e}] vs "
        f"{MAT10_P_PRINT} (band {P_BAND:.0e}); (B) = {err_B:.3e} (band {B_BAND:.0e}); "
        f"(C) = {err_C:.3e} (band {C_BAND:.0e}); RG3:149 {RG3_R}")
    _p0(f"[WF-17 s13] (c-2G) control source_scale {c2.source_scale!r} (2G {2.0 * G!r}): (A) "
        f"comparator {c_2g:.4e} (ASSERTED > {A_BAND}; RG3:150 {RG3_C1}); {c_2g / A_BAND:.1f}x "
        f"the band; no re-solve {once_c}")
    _p0(f"[WF-17 s13] W-heat block (run_study_file) {t_block:.1f} s; node "
        f"{time.perf_counter() - t0:.1f} s")

    assert n_calls == 1
    assert paths == [THERMAL_FILE.name]
    assert s["family"] == "sphere_in_box" and s["degree"] == 1
    assert list(s["temperature"]) == KEYS
    assert list(s["temperature"]) == extract_tagged_block(
        _text(), KEYS_TAG, lang="text").strip("\n").splitlines()
    assert n_sub == MAT10_SUB_CELLS, f"(H) sub cells {n_sub}"
    assert f"{err_A:.4e}" == MAT10_A_PRINT, f"(H) (A) {err_A:.4e}"
    assert f"{err_P:.4e}" == MAT10_P_PRINT, f"(H) (P) {err_P:.4e}"
    assert err_A <= A_BAND and err_P <= P_BAND
    assert err_B <= B_BAND, f"(B) {err_B:.3e}"
    assert err_C <= C_BAND, f"(C) {err_C:.3e}"
    assert once_c, "(c-2G) re-solved"
    assert c_2g > A_BAND, f"(c-2G) {c_2g:.4e} did not miss {A_BAND}"


def test_wf17_s13_probes_block_reproduces_mat19(monkeypatch):
    comm = _solve_preamble("probes")
    t0 = time.perf_counter()
    ns, calls, _ = exec_sphere_block(PROBES_TAG, monkeypatch)
    t_block = time.perf_counter() - t0
    run = ns["run"]
    s = run.summary
    n_calls, paths = calls["n"], list(calls["paths"])
    pr = run.probes
    assert pr is not None
    assert run.study.result is run.result, "(Once) a second EM solve happened"

    # (Pr) on run.probes' θ, as test_wf1_yaml_probes.py does (never run.thermal's).
    sub, theta = pr.thermal.submesh, pr.theta
    cells = int(sub.topology.index_map(3).size_global)
    v_sub = float(np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0))
                                  * ufl.dx(domain=sub))))
    mean = float(np.real(_allsum(comm, theta * ufl.dx(domain=sub)))) / v_sub
    a_val = abs(s["probes"]["values_k"][0] - mean) / abs(mean)
    a_vol = abs(s["probes"]["volumes_m3"][0] / v_sub - 1.0)
    direct = [float(ball_average(pr.theta, p)[0]) for p in pr.probes]
    equal = s["probes"]["values_k"] == direct
    mem_equal = s["probes"]["values_k"] == [float(v) for v in pr.values]

    # (c-c2) c2 displaced 0.1 mm outward misses the summary's c2.
    v_out = float(ball_average(pr.theta, C2_OUT)[0])
    c2v = s["probes"]["values_k"][3]
    c_c2 = abs(v_out - c2v) / abs(c2v)

    _p0(f"[WF-17 s13] (Pr-run) block spy calls {n_calls} {paths}; family {s['family']!r}; "
        f"degree {s['degree']!r}; summary['probes'] {s['probes']}")
    _p0(f"[WF-17 s13] (Pr) sub cells {cells} (N_SUB {N_SUB}); |P - int/V|/|int/V| = {a_val:.3e}, "
        f"|V_h/V_sub - 1| = {a_vol:.3e} (band {B_BAND:.0e}; RGP:153 {RGP_A}); V_sub "
        f"{v_sub:.10e}; whole-ball mean {mean:.10e} K")
    _p0(f"[WF-17 s13] (Pr) summary values {s['probes']['values_k']!r}; ball_average direct "
        f"{direct!r}; equal {equal}; == run.probes.values {mem_equal}")
    _p0(f"[WF-17 s13] (c-c2) c2 +0.1 mm: |d|/theta = {c_c2:.4e} (ASSERTED > {WIDTH_BAND:.0e}; "
        f"RGP:235 {RGP_C1})")
    _p0(f"[WF-17 s13] W-probes block (run_study_file) {t_block:.1f} s; node "
        f"{time.perf_counter() - t0:.1f} s")

    assert n_calls == 1
    assert paths == [PROBES_FILE.name]
    assert s["family"] == "sphere_in_box" and s["degree"] == 1
    assert cells == N_SUB, f"(Pr) {cells}"
    assert a_val <= B_BAND, f"(Pr) value {a_val:.3e}"
    assert a_vol <= B_BAND, f"(Pr) volume {a_vol:.3e}"
    assert equal, (s["probes"]["values_k"], direct)
    assert mem_equal
    assert c_c2 > WIDTH_BAND, f"(c-c2) {c_c2:.3e} did not miss {WIDTH_BAND:.0e}"
