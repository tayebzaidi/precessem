"""`WF-17` step 37 (§9 item 927, §10 U145): the sphere page's default-degree section is a test.

``docs/guide/sphere_heating_transient.md`` gains "At the `Study`'s default
degree": one fenced ``python sphere-transient-degree2`` block (U143's pulsed
block with the ``degree=`` keyword removed, plus a ``study degree`` print),
extracted verbatim and exec'd by item 257's helpers (imported, not copied).
Brief: docs/planning/chunks/WF-17.md, heading "§9 item 927".  Bands imported,
none set; no degree-2 digit is pinned.

W0 (no solve):
(X)  exactly one tagged block, one ``Study(`` with **no** ``degree`` keyword,
     one ``DutyCycle(``, no ``probe``;
(K)  ``Study.solve`` stubbed, ``Study.temperature`` spied (``duty`` recorded from
     ``**kw``): the Study equals ``_study()`` (the default, degree 2) on
     ``STUDY_FIELDS``; cooling ``RobinCooling(K, H, 0.0, source_scale=G)``; both
     specs ``EXPECTED_A``; call 1 no duty, call 2 ``EXPECTED_DUTY``;
(E)  the ``sphere-transient-degree2-output`` block is exactly
     ``sphere cells = MAT10_SUB_CELLS`` and ``study degree = 2``;
(D)  docrefs dead 0 and ``MIN_GUIDE_PAGES`` holds.
W1 (``FEM_EM_WF17_S37=1``): (Deg), (Once), (UT), (χ), (S0), (Δ) ≤ ``L_BAND``,
(A_t) ≤ ``A_BAND`` (``.4e`` printed beside U132's degree-2 record and MAT-12's
degree-1 record, never asserted equal), (E) ≤ ``E_BAND`` on 12 runs.
Controls (asserted): (Δ-ctl) dropping m = 9 misses ``L_BAND``; (c-sched) the
pulsed θ misses ``A_BAND`` against the continuous series; (c-deg) a planted
``, degree=1`` makes (K)'s mismatch list exactly ``['study.degree']`` and fires
(X)'s keyword check; (c-page) a removed or doubled tag raises ``LookupError``.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pytest
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver import DutyCycle, RobinCooling, Study
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.package.test_wf15_study_sphere_transient import (
    A_BAND,
    E_BAND,
    G,
    H,
    K,
    L_BAND,
    MAT10_SUB_CELLS,
    MAT12_A_T_PRINT,
    N_MAIN,
    T_STAR,
    U132_A_T_DEG2_PRINT,
    U135_CHI_PRINT,
    _census,
    _chi_string,
    _exact_theta,
    _on_space,
    _space_degree,
    _study,
    _whole_rel,
)
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    _study_calls,
    exec_block,
    extract_tagged_block,
)
from tests.package.test_wf17_sphere_transient_duty_tutorial import EXPECTED_DUTY
from tests.package.test_wf17_sphere_transient_tutorial import (
    COOL_FIELDS,
    EXPECTED_A,
    MIN_GUIDE_PAGES,
    PAGE,
    SPEC_FIELDS,
    STUDY_FIELDS,
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
    ROBIN_TAG,
    _interior_indicator,
    _interior_rel,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C

TAG = "sphere-transient-degree2"
OUTPUT_TAG = "sphere-transient-degree2-output"
REPO_ROOT = Path(__file__).resolve().parents[2]
SCRATCH = REPO_ROOT / "logs" / "wf17" / "s37"
S37_ENV = "FEM_EM_WF17_S37"
# MAT-12's degree-1 (A_t) record, printed beside the degree-2 reading (c-deg).
DEG1_A_T_PRINT = MAT12_A_T_PRINT

needs_s37 = pytest.mark.skipif(
    os.environ.get(S37_ENV) != "1", reason=f"W1 solve: set {S37_ENV}=1"
)


def _say(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _page_and_block():
    text = PAGE.read_text()
    return text, extract_tagged_block(text, TAG)


def _degree_keywords(block) -> list[str]:
    """(X)'s keyword check: every ``degree`` keyword on the block's ``Study(`` calls."""
    return [kw.arg for call in _study_calls(block) for kw in call.keywords
            if kw.arg == "degree"]


class _StopAfterSpecs(Exception):
    pass


def _mismatches(study, calls) -> list[str]:
    """(K)'s comparator against the default ``_study()`` (degree 2)."""
    ref = _study()
    out = [f"study.{n}" for n in STUDY_FIELDS if getattr(study, n) != getattr(ref, n)]
    if len(calls) != 2:
        return out + [f"{len(calls)} temperature calls"]
    cool = RobinCooling(K, H, 0.0, source_scale=G)
    for label, (cooling, spec, duty), exp_duty in zip("AD", calls, (None, EXPECTED_DUTY)):
        out += [f"{label}.cooling.{n}" for n in COOL_FIELDS
                if getattr(cooling, n) != getattr(cool, n)]
        if spec is None:
            out.append(f"{label}.transient None")
        else:
            out += [f"{label}.{n}" for n in SPEC_FIELDS
                    if getattr(spec, n) != getattr(EXPECTED_A, n)]
        if exp_duty is None:
            if duty is not None:
                out.append(f"{label}.duty not None")
        elif not isinstance(duty, DutyCycle):
            out.append(f"{label}.duty {duty!r}")
        else:
            out += [f"{label}.duty.{n}" for n in ("period_s", "duty")
                    if getattr(duty, n) != getattr(exp_duty, n)]
    return out


def _exec_no_solve(block, monkeypatch):
    """Exec with Study.solve stubbed and Study.temperature spied (records ``duty``)."""
    rec: dict = {"calls": []}

    def stub_solve(self, comm=MPI.COMM_WORLD):
        rec["study"] = self
        return object()

    def spy_temperature(self, cooling, transient=None, **kw):
        rec["study"] = self
        rec["calls"].append((cooling, transient, kw.get("duty")))
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


def test_wf17_s37_extracts_exactly_one_block():
    _, block = _page_and_block()
    calls = _study_calls(block)
    n_duty = block.count("DutyCycle(")
    deg_kw = _degree_keywords(block)
    _say(f"[WF-17s37] (X) Study( calls {len(calls)}; degree keywords {deg_kw}; "
         f"DutyCycle( {n_duty}; probe in block: {'probe' in block}")
    assert len(calls) == 1
    assert deg_kw == []
    assert n_duty == 1
    assert "probe" not in block


def test_wf17_s37_block_configures_the_default(monkeypatch):
    _, block = _page_and_block()
    rec = _exec_no_solve(block, monkeypatch)
    mism = _mismatches(rec["study"], rec["calls"])
    duties = [None if d is None else f"DutyCycle({d.period_s:.6g}, {d.duty})"
              for d in (c[2] for c in rec["calls"])]
    _say(f"[WF-17s37] (K) degree={rec['study'].degree} (_study() degree "
         f"{_study().degree}) duties={duties} mismatches={mism}")
    assert rec["study"].degree == 2
    assert mism == []


def test_wf17_s37_planted_degree1_is_caught(monkeypatch):
    _, block = _page_and_block()
    old = "density_kg_m3=1000.0)"
    assert block.count(old) == 1
    planted = block.replace(old, "density_kg_m3=1000.0, degree=1)")
    deg_kw = _degree_keywords(planted)
    rec = _exec_no_solve(planted, monkeypatch)
    mism = _mismatches(rec["study"], rec["calls"])
    _say(f"[WF-17s37] (c-deg) planted ', degree=1': mismatches={mism}; "
         f"(X) degree keywords {deg_kw}")
    assert mism == ["study.degree"]
    assert deg_kw == ["degree"]


def test_wf17_s37_planted_untagged_and_doubled_pages_raise():
    text, block = _page_and_block()
    untagged = text.replace(f"```python {TAG}\n", "```python\n")
    with pytest.raises(LookupError):
        extract_tagged_block(untagged, TAG)
    doubled = text + f"\n```python {TAG}\n{block}```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    _say("[WF-17s37] (c-page) untagged and doubled pages both raise LookupError")


def test_wf17_s37_expected_output_is_two_lines():
    text, _ = _page_and_block()
    out = extract_tagged_block(text, OUTPUT_TAG, lang="text")
    expected = f"sphere cells = {MAT10_SUB_CELLS}\nstudy degree = 2\n"
    _say(f"[WF-17s37] (E) {out!r} vs {expected!r}")
    assert out == expected


def test_wf17_s37_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17s37] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"(MIN_GUIDE_PAGES {MIN_GUIDE_PAGES}) exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ---------------------------------------------------------------- W1 (the solve)


@complex_only
@needs_s37
def test_wf17_s37_block_at_default_degree(monkeypatch):
    comm = MPI.COMM_WORLD
    _, block = _page_and_block()
    rec: dict = {"calls": [], "thermal": []}
    orig_solve, orig_temp = Study.solve, Study.temperature

    def spy_solve(self, comm=MPI.COMM_WORLD):
        rec["study"] = self
        t = time.perf_counter()
        rec["result"] = orig_solve(self, comm)
        rec["t_em"] = time.perf_counter() - t
        return rec["result"]

    def spy_temperature(self, cooling, transient=None, **kw):
        th = orig_temp(self, cooling, transient=transient, **kw)
        rec["calls"].append((cooling, transient, kw.get("duty")))
        rec["thermal"].append(th)
        return th

    monkeypatch.setattr(Study, "solve", spy_solve)
    monkeypatch.setattr(Study, "temperature", spy_temperature)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(SCRATCH)
    t0 = time.perf_counter()
    _, captured = exec_block(block, monkeypatch, stop_before_solve=False)
    t_block = time.perf_counter() - t0

    study, res = rec["study"], rec["result"]
    mism = _mismatches(study, rec["calls"])
    heating, pulsed = rec["thermal"]
    once = study.result is res

    # (Deg)
    sp_deg = _space_degree(res.fields.e_complex.function_space)
    _say(f"WF17S37 (Deg) study.degree {study.degree}, res.degree {res.degree}, "
         f"E space degree {sp_deg} (U144 read 2 / 2 / 2); EM solve {rec['t_em']:.1f} s")

    dt = T_STAR / N_MAIN
    t1 = time.perf_counter()
    H_runs = []
    for m in range(10):
        H_runs.append(solve_transient_conduction(
            pulsed.mesh, pulsed.facet_tags, k=K, rho_c=RHO_C, q=pulsed.q, h=H, t_inf=0.0,
            t_init=0.0, robin_tag=ROBIN_TAG, dt=dt, n_steps=N_MAIN, theta=0.5,
            startup_half_steps=2, degree=1,
            source_schedule=lambda t, tm=m * 4 * dt: 1.0 if t > tm else 0.0,
        ))
    t_dir = time.perf_counter() - t1
    once = once and study.result is res

    # (UT)
    n_sub = int(heating.mesh.topology.index_map(3).size_global)
    n_sub_p = int(pulsed.mesh.topology.index_map(3).size_global)
    _say(f"WF17S37 (K) mismatches={mism}; once {once}; sub cells heating {n_sub} "
         f"pulsed {n_sub_p} (MAT10_SUB_CELLS {MAT10_SUB_CELLS})")

    # (χ) / (S0)
    chis = list(pulsed.source_multipliers)
    chi_binary = all(c in (0.0, 1.0) for c in chis)
    chi_s = _chi_string(chis)
    _say(f"WF17S37 (chi) chi string {chi_s} (n {len(chis)}; binary {chi_binary}; "
         f"== U135_CHI_PRINT {chi_s == U135_CHI_PRINT})")
    _census(pulsed, dt, 2, T_STAR / 2, "WF17S37 pulsed")

    # (Δ) Duhamel on pulsed's space.
    sub = pulsed.mesh
    V1 = pulsed.temperature.function_space
    lay_H = True
    S = fem.Function(V1)
    S_ctl = fem.Function(V1)
    for m, h_run in enumerate(H_runs):
        arr = h_run.temperature.x.array
        ok = bool(comm.allreduce(bool(arr.size == S.x.array.size), op=MPI.LAND))
        lay_H = lay_H and ok
        if ok:
            S.x.array[:] += ((-1) ** m) * arr
            if m != 9:
                S_ctl.x.array[:] += ((-1) ** m) * arr
    err_L = _whole_rel(sub, pulsed.temperature - S, pulsed.temperature)
    err_L_ctl = _whole_rel(sub, pulsed.temperature - S_ctl, pulsed.temperature)
    _say(f"WF17S37 (Delta) whole rel L2 theta_pulsed - sum(-1)^m theta_H_m = {err_L:.3e} "
         f"(band {L_BAND:.0e}; U144 read 6.638e-16); layouts {lay_H}")
    _say(f"WF17S37 (Delta-ctl) dropping m = 9: {err_L_ctl:.6f} (must miss {L_BAND:.0e}; "
         f"U144 read 0.284784)")

    # (A_t) on heating's own submesh; (c-sched) on pulsed's.
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    subA = heating.mesh
    theta_A = _exact_theta(subA, T_STAR, q_cf)
    err_A = _interior_rel(subA, _interior_indicator(subA), heating.temperature - theta_A,
                          theta_A)
    _say(f"WF17S37 (A_t) heating interior rel L2 = {err_A:.10e} [{err_A:.4e}] "
         f"(band {A_BAND}; printed, not pinned: U132 degree-2 record {U132_A_T_DEG2_PRINT}; "
         f"degree-1 record {DEG1_A_T_PRINT}; ratio deg1/deg2 "
         f"{float(DEG1_A_T_PRINT) / err_A:.1f})")
    theta_cf = _exact_theta(sub, T_STAR, q_cf)
    err_cs = _interior_rel(sub, _interior_indicator(sub), pulsed.temperature - theta_cf,
                           theta_cf)
    _say(f"WF17S37 (c-sched) theta_pulsed vs continuous series: interior rel L2 = "
         f"{err_cs:.4e} (must miss {A_BAND}; U144 read 5.1677e-01)")

    _, lay_A = _on_space(V1, heating.temperature, comm)

    # (E)
    runs = [("heating", heating), ("pulsed", pulsed)] + [
        (f"H{m}", h) for m, h in enumerate(H_runs)]
    es = [(lbl, float(r.max_balance_residual)) for lbl, r in runs]
    _say("WF17S37 (E) max balance " + "  ".join(f"{lbl} {e:.3e}" for lbl, e in es)
         + f"  (band {E_BAND:.0e}; {len(es)} runs)")
    lines = [f"sphere cells = {MAT10_SUB_CELLS}", "study degree = 2"]
    _say(f"WF17S37 heating/pulsed same layout+geometry {lay_A}; stdout has "
         + ", ".join(f"{ln!r}: {ln in captured}" for ln in lines))
    _say(f"WF17S37 elapsed: block (1 EM + 2 transients) {t_block:.1f} s, direct runs (10) "
         f"{t_dir:.1f} s, total {time.perf_counter() - t0:.1f} s")

    assert study.degree == 2 and res.degree == 2 and sp_deg == 2, \
        f"(Deg) {study.degree} {res.degree} {sp_deg}"
    assert mism == [], f"(K) {mism}"
    assert once, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS and n_sub_p == MAT10_SUB_CELLS, f"(UT) {n_sub} {n_sub_p}"
    assert chi_binary, f"(chi) {chi_s}"
    assert chi_s == U135_CHI_PRINT, f"(chi) {chi_s}"
    assert lay_H, "(Delta) layouts differ"
    assert err_L <= L_BAND, f"(Delta) {err_L:.3e} > {L_BAND:.0e}"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > {A_BAND}"
    assert all(e <= E_BAND for _, e in es), f"(E) {es}"
    assert err_L_ctl > L_BAND, f"(Delta-ctl) {err_L_ctl:.3e} did not miss {L_BAND:.0e}"
    assert err_cs > A_BAND, f"(c-sched) {err_cs:.4e} did not miss {A_BAND}"
    if comm.rank == 0:
        assert all(ln in captured for ln in lines), "(E) output lines missing"
    print(f"WF17S37 VERDICT rank {comm.rank} PASS", flush=True)
