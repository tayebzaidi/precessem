"""`WF-17` step 33 (§9 item 865, §10 U133): the SAR-driven sphere transient is a test.

``docs/guide/sphere_heating_transient.md`` (tag ``sphere-transient``, from
`EX-187` / ``examples/mri/40_study_sar_driven_transient_sphere.py``) carries one
fenced ``python sphere-transient`` block, extracted verbatim and exec'd by item
257's helpers (imported from ``test_wf17_first_study_tutorial``, not copied).
Asserted (brief: docs/planning/chunks/WF-17.md, heading "§9 item 865"):

(X)  exactly one tagged block, with exactly one ``Study(`` call;
(K)  no solve: the block's Study equals ``_study(degree=1)`` field by field, its
     ``RobinCooling`` is ``(K, H, 0.0, source_scale=G)``, and its two
     ``TransientSpec``s are U131's A and S;
(E)  the ``sphere-transient-output`` block states ``MAT12_A_T_PRINT``,
     ``MAT12_S2_PRINT`` and ``MAT10_SUB_CELLS`` literally;
(R)  W1 (``FEM_EM_WF17_S33=1``): (A_t) <= ``A_BAND`` and ``.4e`` ==
     ``MAT12_A_T_PRINT``; (S2) <= ``S2_BAND`` and ``.4e`` == ``MAT12_S2_PRINT``;
     sub cells == ``MAT10_SUB_CELLS``; (E) <= ``E_BAND``; ``study.result is res``;
(D)  ``docrefs --docs-root docs/guide --check-links`` reads dead == 0 and
     ``MIN_GUIDE_PAGES`` holds.

Controls (solve-free, asserted): (c-page) removed or doubled tag raises
``LookupError``; (c-K) the planted edit ``n_steps = 40 -> 20`` makes (K)'s
mismatch list exactly the A spec's ``dt_s`` and ``n_steps``.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import RobinCooling, Study, TransientSpec
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import _study
from tests.package.test_wf15_study_sphere_transient import (
    MAT12_A_T_PRINT,
    MAT12_S2_PRINT,
)
from tests.package.test_wf15_study_temperature import MAT10_SUB_CELLS
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    _study_calls,
    exec_block,
    extract_tagged_block,
)
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
    G,
    H,
    K,
    _allsum,
    _interior_indicator,
    _interior_rel,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    DT_1,
    E_BAND,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_STAR,
    _exact_theta,
    _mean,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
PAGE = REPO_ROOT / "docs" / "guide" / "sphere_heating_transient.md"
TAG = "sphere-transient"
OUTPUT_TAG = "sphere-transient-output"
SCRATCH = REPO_ROOT / "logs" / "wf17" / "s33"
S33_ENV = "FEM_EM_WF17_S33"
# W0 of this step read 76 guide pages with the new page (item 838 read 75 before it).
MIN_GUIDE_PAGES = 76

STUDY_FIELDS = ("geometry", "material", "excitation", "frequency_hz", "degree",
                "density_kg_m3", "conductor")
SPEC_FIELDS = ("dt_s", "n_steps", "rho_c_j_per_m3_k", "theta", "startup_half_steps",
               "t_init_c", "initial_gradient_k_per_m")
COOL_FIELDS = ("k_w_per_m_k", "h_w_per_m2_k", "t_inf_c", "source_scale", "degree")
EXPECTED_A = TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)
EXPECTED_S = TransientSpec(DT_1, 1, RHO_C, theta=1.0, startup_half_steps=0)

needs_s33 = pytest.mark.skipif(
    os.environ.get(S33_ENV) != "1", reason=f"W1 solve: set {S33_ENV}=1"
)


def _say(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _page_and_block():
    text = PAGE.read_text()
    return text, extract_tagged_block(text, TAG)


class _StopAfterSpecs(Exception):
    pass


def _mismatches(study, calls) -> list[str]:
    """(K)'s comparator: Study fields, the cooling, and the two TransientSpecs."""
    ref = _study(degree=1)
    out = [f"study.{n}" for n in STUDY_FIELDS if getattr(study, n) != getattr(ref, n)]
    if len(calls) != 2:
        return out + [f"{len(calls)} temperature calls"]
    cool = RobinCooling(K, H, 0.0, source_scale=G)
    for label, (cooling, spec), exp in zip("AS", calls, (EXPECTED_A, EXPECTED_S)):
        out += [f"{label}.cooling.{n}" for n in COOL_FIELDS
                if getattr(cooling, n) != getattr(cool, n)]
        if spec is None:
            out.append(f"{label}.transient None")
            continue
        out += [f"{label}.{n}" for n in SPEC_FIELDS if getattr(spec, n) != getattr(exp, n)]
    return out


def _exec_no_solve(block, monkeypatch):
    """Exec with Study.solve stubbed and Study.temperature spied (no solve at all)."""
    rec: dict = {"calls": []}

    def stub_solve(self, comm=MPI.COMM_WORLD):
        rec["study"] = self
        return object()

    def spy_temperature(self, cooling, transient=None, **kw):
        rec["study"] = self
        rec["calls"].append((cooling, transient))
        if len(rec["calls"]) == 2:
            raise _StopAfterSpecs
        return object()

    monkeypatch.setattr(Study, "solve", stub_solve)
    monkeypatch.setattr(Study, "temperature", spy_temperature)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(SCRATCH)
    with pytest.raises(_StopAfterSpecs):
        exec_block(block, monkeypatch, stop_before_solve=False)
    return rec


# ---------------------------------------------------------------- W0 (no solve)


def test_wf17_s33_extracts_exactly_one_block():
    _, block = _page_and_block()
    calls = _study_calls(block)
    keywords = {kw.arg: kw.value for kw in calls[0].keywords} if calls else {}
    _say(f"[WF-17s33] (X) Study( calls {len(calls)}; keywords {sorted(keywords)}; "
         f"probe in block: {'probe' in block}")
    assert len(calls) == 1
    assert "degree" in keywords
    assert "probe" not in block


def test_wf17_s33_block_configures_u131(monkeypatch):
    _, block = _page_and_block()
    rec = _exec_no_solve(block, monkeypatch)
    mism = _mismatches(rec["study"], rec["calls"])
    _say(f"[WF-17s33] (K) degree={rec['study'].degree} specs={[c[1] for c in rec['calls']]} "
         f"mismatches={mism}")
    assert rec["study"].degree == 1
    assert mism == []


def test_wf17_s33_planted_n_steps_is_caught(monkeypatch):
    _, block = _page_and_block()
    assert block.count("n_steps = 40\n") == 1
    planted = block.replace("n_steps = 40\n", f"n_steps = {N_MAIN // 2}\n")
    rec = _exec_no_solve(planted, monkeypatch)
    mism = _mismatches(rec["study"], rec["calls"])
    _say(f"[WF-17s33] (c-K) n_steps -> {N_MAIN // 2}: mismatches={mism}")
    assert sorted(mism) == ["A.dt_s", "A.n_steps"]


def test_wf17_s33_planted_untagged_and_doubled_pages_raise():
    text, block = _page_and_block()
    untagged = text.replace(f"```python {TAG}\n", "```python\n")
    with pytest.raises(LookupError):
        extract_tagged_block(untagged, TAG)
    doubled = text + f"\n```python {TAG}\n{block}```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    _say("[WF-17s33] (c-page) untagged and doubled pages both raise LookupError")


def test_wf17_s33_expected_output_states_the_records():
    text, _ = _page_and_block()
    out = extract_tagged_block(text, OUTPUT_TAG, lang="text")
    expected = (f"sphere cells = {MAT10_SUB_CELLS}\n(A_t) {MAT12_A_T_PRINT}\n"
                f"(S2) {MAT12_S2_PRINT}\n")
    _say(f"[WF-17s33] (E) {out!r} vs {expected!r}")
    assert out == expected


def test_wf17_s33_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17s33] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"(MIN_GUIDE_PAGES {MIN_GUIDE_PAGES}) exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ---------------------------------------------------------------- W1 (the solve)


@complex_only
@needs_s33
def test_wf17_s33_block_reproduces_u131(monkeypatch):
    comm = MPI.COMM_WORLD
    _, block = _page_and_block()
    rec: dict = {"calls": []}
    orig_solve, orig_temp = Study.solve, Study.temperature

    def spy_solve(self, comm=MPI.COMM_WORLD):
        rec["study"] = self
        rec["result"] = orig_solve(self, comm)
        return rec["result"]

    def spy_temperature(self, cooling, transient=None, **kw):
        th = orig_temp(self, cooling, transient=transient, **kw)
        rec["calls"].append((cooling, transient))
        rec.setdefault("thermal", []).append(th)
        return th

    monkeypatch.setattr(Study, "solve", spy_solve)
    monkeypatch.setattr(Study, "temperature", spy_temperature)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(SCRATCH)
    t0 = time.perf_counter()
    _, captured = exec_block(block, monkeypatch, stop_before_solve=False)
    seconds = time.perf_counter() - t0

    study, res = rec["study"], rec["result"]
    mism = _mismatches(study, rec["calls"])
    thA, thS = rec["thermal"]
    once = study.result is res

    n_sub = int(thA.mesh.topology.index_map(3).size_global)
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub = thA.mesh
    theta = _exact_theta(sub, T_STAR, q_cf)
    err_A = _interior_rel(sub, _interior_indicator(sub), thA.temperature - theta, theta)

    subS = thS.mesh
    one = fem.Constant(subS, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=subS))))
    ratio = RHO_C * _mean(subS, thS.temperature, vol) / (DT_1 * _mean(subS, thS.q, vol))
    s2 = abs(ratio - 1.0)
    eA, eS = thA.max_balance_residual, thS.max_balance_residual

    _say(f"WF17S33 (K) mismatches={mism}; once {once}; sub cells {n_sub} "
         f"(MAT10_SUB_CELLS {MAT10_SUB_CELLS})")
    _say(f"WF17S33 (A_t) {err_A:.10e} [{err_A:.4e}] (band {A_BAND}; U131 {MAT12_A_T_PRINT})")
    _say(f"WF17S33 (S2) ratio {ratio:.9f} |ratio - 1| {s2:.4e} (band {S2_BAND:.0e}; "
         f"U131 {MAT12_S2_PRINT})")
    _say(f"WF17S33 (E) max balance A {eA:.3e} S {eS:.3e} (band {E_BAND:.0e})")
    _say(f"WF17S33 block exec {seconds:.1f} s")

    assert mism == [], f"(K) {mism}"
    assert once, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(UT) sub cells {n_sub}"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > {A_BAND}"
    assert f"{err_A:.4e}" == MAT12_A_T_PRINT, f"(A_t) print {err_A:.4e}"
    assert s2 <= S2_BAND, f"(S2) {s2:.4e} > {S2_BAND:.0e}"
    assert f"{s2:.4e}" == MAT12_S2_PRINT, f"(S2) print {s2:.4e}"
    assert eA <= E_BAND and eS <= E_BAND, f"(E) A {eA:.3e} S {eS:.3e}"
    if comm.rank == 0:
        assert f"sphere cells = {MAT10_SUB_CELLS}" in captured
