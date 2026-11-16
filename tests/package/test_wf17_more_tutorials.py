"""`WF-17` step 3 (§9 item 329, chain U29): two more tutorials are tests.

``docs/guide/four_port_s.md`` (tag ``four-port-s``, from `EX-92`) and
``docs/guide/temperature.md`` (tag ``temperature``, from `EX-94`) each carry one
fenced ``python <tag>`` block, extracted verbatim and exec'd by item 257's
helpers (imported from ``test_wf17_first_study_tutorial``, not copied). Asserted
(brief: docs/planning/chunks/WF-17.md, heading "§9 item 329"):

(X)  exactly one tagged block per page, with exactly one ``Study(`` call;
(K)  each block configures the gate's Study exactly, at ``degree == 1``
     (no solve: ``Study.solve`` is stubbed);
(E)  each page's ``text <tag>-output`` block states the record literally;
(D)  ``docrefs --docs-root docs/guide --check-links`` reads dead == 0;
(R4) four-port: ``_rel(S, P1) <= LEG_D0_REPRODUCTION_BAND`` at the record width;
(A)/(P) temperature: equal `MAT-10`'s records at ``.4e``; (A)/(P)/(B)/(C) inside
     their imported bands; sphere cells == ``MAT10_SUB_CELLS``.

Negative controls (solve-free, asserted): (c1) untagged / doubly-tagged page
raises, per page; (c2) a planted parameter edit (frequency 64 -> 128 MHz; h ->
10 h) makes (K)'s list exactly the edited field; (c3) a planted dead link in
``tmp_path`` reads dead == 1 (``test_wf17_docrefs_detects_a_planted_dead_link``
in the step-1 module, run in W0).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver import (
    FourPortBirdcageGeometry,
    GelledSalinePhantomMaterial,
    LumpedPortDrive,
    QuasiStaticSphereField,
    RobinCooling,
    Study,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import GEOMETRY
from tests.package.test_wf15_study_temperature import (
    MAT10_A_PRINT,
    MAT10_P_PRINT,
    MAT10_SUB_CELLS,
)
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    _study_calls,
    exec_block,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    FREQUENCY_HZ,
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
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

REPO_ROOT = Path(__file__).resolve().parents[2]
GUIDE = REPO_ROOT / "docs" / "guide"
FOUR_PORT_PAGE = GUIDE / "four_port_s.md"
TEMPERATURE_PAGE = GUIDE / "temperature.md"
FOUR_PORT_TAG = "four-port-s"
TEMPERATURE_TAG = "temperature"
PAGES = {FOUR_PORT_TAG: FOUR_PORT_PAGE, TEMPERATURE_TAG: TEMPERATURE_PAGE}


def _say(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _page_and_block(tag):
    text = PAGES[tag].read_text()
    return text, extract_tagged_block(text, tag)


def _assert_extraction(block):
    assert block.strip(), "(X) empty block"
    calls = _study_calls(block)
    assert len(calls) == 1, f"(X) {len(calls)} Study( calls in the block, need exactly 1"
    return calls[0]


def _four_port_mismatches(study) -> list[str]:
    expected = {
        "geometry": FourPortBirdcageGeometry.f_small(),
        "material": GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ),
        "excitation": LumpedPortDrive(),
        "frequency_hz": FREQUENCY_64_HZ,
        "degree": 1,
        "conductor": None,
    }
    return [name for name, value in expected.items() if getattr(study, name) != value]


def _temperature_mismatches(study, cooling) -> list[str]:
    expected = {
        "geometry": GEOMETRY,
        "material": GelledSalinePhantomMaterial(SIGMA_HIGH, EPSILON_R_SPHERE, FREQUENCY_HZ),
        "excitation": QuasiStaticSphereField(e0=E0),
        "frequency_hz": FREQUENCY_HZ,
        "density_kg_m3": RHO_KG_M3,
        "degree": 1,
    }
    out = [name for name, value in expected.items() if getattr(study, name) != value]
    spec = RobinCooling(K, H, T_INF, source_scale=G)
    for name in ("k_w_per_m_k", "h_w_per_m2_k", "t_inf_c", "source_scale", "degree"):
        if getattr(cooling, name) != getattr(spec, name):
            out.append(f"cooling.{name}")
    return out


def _spy_temperature(monkeypatch, recorded):
    original = Study.temperature

    def spy(self, cooling):
        recorded["study"] = self
        recorded["cooling"] = cooling
        th = original(self, cooling)
        recorded["thermal"] = th
        return th

    monkeypatch.setattr(Study, "temperature", spy)


def _spy_solve(monkeypatch, recorded):
    original = Study.solve

    def spy(self, comm=MPI.COMM_WORLD):
        recorded["study"] = self
        res = original(self, comm)
        recorded["result"] = res
        return res

    monkeypatch.setattr(Study, "solve", spy)


def _exec_temperature(block, monkeypatch, *, stop_before_solve):
    """Exec the temperature block; returns (spy dict, exec's dict, stdout)."""
    spied: dict = {}
    _spy_temperature(monkeypatch, spied)
    recorded, captured = exec_block(block, monkeypatch, stop_before_solve=stop_before_solve)
    return spied, recorded, captured


# ---------------------------------------------------------------- W0 (smoke)


@pytest.mark.parametrize("tag", [FOUR_PORT_TAG, TEMPERATURE_TAG])
def test_wf17_step3_extracts_exactly_one_block(tag):
    _, block = _page_and_block(tag)
    call = _assert_extraction(block)
    keywords = {kw.arg: kw.value for kw in call.keywords}
    _say(f"[WF-17s3] {tag}: Study keywords {sorted(keywords)}")
    assert "degree" in keywords, "(K) the block must set degree explicitly (degree trap)"


def test_wf17_four_port_block_configures_the_p1_study(monkeypatch):
    _, block = _page_and_block(FOUR_PORT_TAG)
    recorded, _ = exec_block(block, monkeypatch, stop_before_solve=True)
    assert "study" in recorded, "(K) the block never reached Study.solve"
    mismatches = _four_port_mismatches(recorded["study"])
    _say(f"[WF-17s3] (K) four-port degree={recorded['study'].degree} mismatches={mismatches}")
    assert recorded["study"].degree == 1
    assert mismatches == []


def test_wf17_temperature_block_configures_the_mat10_study(monkeypatch):
    _, block = _page_and_block(TEMPERATURE_TAG)
    spied, recorded, _ = _exec_temperature(block, monkeypatch, stop_before_solve=True)
    assert "cooling" in spied and "study" in recorded, "(K) the block never reached the solve"
    assert "thermal" not in spied
    mismatches = _temperature_mismatches(recorded["study"], spied["cooling"])
    _say(f"[WF-17s3] (K) temperature degree={recorded['study'].degree} mismatches={mismatches}")
    assert recorded["study"].degree == 1
    assert mismatches == []


@pytest.mark.parametrize("tag", [FOUR_PORT_TAG, TEMPERATURE_TAG])
def test_wf17_step3_planted_untagged_and_doubled_pages_raise(tag):
    text, block = _page_and_block(tag)
    untagged = text.replace(f"```python {tag}\n", "```python\n")
    with pytest.raises(LookupError):
        extract_tagged_block(untagged, tag)
    doubled = text + f"\n```python {tag}\n{block}```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, tag)
    _say(f"[WF-17s3] (c1) {tag}: untagged and doubled pages both raise LookupError")


def test_wf17_four_port_planted_frequency_is_caught(monkeypatch):
    _, block = _page_and_block(FOUR_PORT_TAG)
    assert block.count("frequency_hz=64e6") == 1
    planted = block.replace("frequency_hz=64e6", "frequency_hz=128e6")
    recorded, _ = exec_block(planted, monkeypatch, stop_before_solve=True)
    mismatches = _four_port_mismatches(recorded["study"])
    _say(f"[WF-17s3] (c2) four-port frequency 128 MHz: mismatches={mismatches}")
    assert mismatches == ["frequency_hz"]


def test_wf17_temperature_planted_h_is_caught(monkeypatch):
    _, block = _page_and_block(TEMPERATURE_TAG)
    assert block.count("100.0") == 1
    planted = block.replace("100.0", "1000.0")
    spied, recorded, _ = _exec_temperature(planted, monkeypatch, stop_before_solve=True)
    mismatches = _temperature_mismatches(recorded["study"], spied["cooling"])
    _say(f"[WF-17s3] (c2) temperature h -> 10 h: mismatches={mismatches}")
    assert mismatches == ["cooling.h_w_per_m2_k"]


def test_wf17_step3_expected_outputs_state_the_records():
    p1 = _record(1)
    four_text, _ = _page_and_block(FOUR_PORT_TAG)
    four_out = extract_tagged_block(four_text, f"{FOUR_PORT_TAG}-output", lang="text")
    four_expected = f"|S11| = {abs(p1[0, 0]):.4f}\n|S21| = {abs(p1[1, 0]):.4f}\n"
    temp_text, _ = _page_and_block(TEMPERATURE_TAG)
    temp_out = extract_tagged_block(temp_text, f"{TEMPERATURE_TAG}-output", lang="text")
    temp_expected = f"sphere cells = {MAT10_SUB_CELLS}\n"
    _say(f"[WF-17s3] (E) four-port {four_out!r} vs {four_expected!r}")
    _say(f"[WF-17s3] (E) temperature {temp_out!r} vs {temp_expected!r}; (A) print "
         f"{MAT10_A_PRINT} in page prose: {MAT10_A_PRINT in temp_text}")
    assert four_out == four_expected
    assert temp_out == temp_expected
    assert MAT10_A_PRINT in temp_text


def test_wf17_step3_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17s3] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= 38


# ------------------------------------------------------- complex: W1 and W2


@complex_only
def test_wf17_four_port_block_reproduces_p1(monkeypatch):
    comm = MPI.COMM_WORLD
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}"
    )
    _, block = _page_and_block(FOUR_PORT_TAG)
    spied: dict = {}
    _spy_solve(monkeypatch, spied)
    t0 = time.perf_counter()
    _, captured = exec_block(block, monkeypatch, stop_before_solve=False)
    seconds = time.perf_counter() - t0
    study, res = spied["study"], spied["result"]
    assert _four_port_mismatches(study) == []
    assert study.degree == 1 and res.degree == 1
    s = np.asarray(res.s_matrix, dtype=np.complex128)
    rel = _rel(s, _record(1))
    _say(f"[WF-17s3] (R4) ||S - P1||/||P1|| = {rel:.3e} (band LEG_D0_REPRODUCTION_BAND "
         f"{LEG_D0_REPRODUCTION_BAND:g}); cells {res.cells}; block exec {seconds:.1f} s")
    assert rel <= LEG_D0_REPRODUCTION_BAND
    if comm.rank == 0:
        assert f"|S11| = {abs(s[0, 0]):.4f}" in captured


@complex_only
def test_wf17_temperature_block_reproduces_mat10(monkeypatch):
    comm = MPI.COMM_WORLD
    _, block = _page_and_block(TEMPERATURE_TAG)
    t0 = time.perf_counter()
    spied, _, captured = _exec_temperature(block, monkeypatch, stop_before_solve=False)
    seconds = time.perf_counter() - t0
    study, th = spied["study"], spied["thermal"]
    assert _temperature_mismatches(study, spied["cooling"]) == []
    assert study.degree == 1 and study.result.degree == 1
    res = study.result
    # (A) / (P) / (B) / (C), term for term as EX-94 (examples/mri/10_study_temperature.py).
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub, ft = th.submesh, th.facet_tags
    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    err_a = _interior_rel(sub, chi, (th.temperature - T_INF) - theta, theta)
    res_p = solve_steady_conduction(
        sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_p = _interior_rel(sub, chi, (res_p.temperature - T_INF) - theta, theta)
    int_q = float(np.real(_allsum(comm, th.heat_source * ufl.dx(domain=sub))))
    p_parent = spied["cooling"].source_scale * study.mean_sar()["dissipated_power_w"]
    err_b = abs(int_q - p_parent) / abs(p_parent)
    err_c = abs(th.convected_power - th.source_power) / abs(th.source_power)
    once = study.result is res
    n_sub = int(sub.topology.index_map(sub.topology.dim).size_global)
    _say(f"[WF-17s3] (A) {err_a:.4e} vs {MAT10_A_PRINT} ({err_a:.10e}, band {A_BAND}); "
         f"(P) {err_p:.4e} vs {MAT10_P_PRINT} (band {P_BAND:.0e}); (B) {err_b:.3e} "
         f"(band {B_BAND:.0e}); (C) {err_c:.3e} (band {C_BAND:.0e}); sub cells {n_sub} "
         f"(record {MAT10_SUB_CELLS}); once {once}; block exec {seconds:.1f} s")
    assert once
    assert n_sub == MAT10_SUB_CELLS
    assert f"{err_a:.4e}" == MAT10_A_PRINT
    assert f"{err_p:.4e}" == MAT10_P_PRINT
    assert err_a <= A_BAND
    assert err_p <= P_BAND
    assert err_b <= B_BAND
    assert err_c <= C_BAND
    if comm.rank == 0:
        assert f"sphere cells = {MAT10_SUB_CELLS}" in captured
