"""`WF-17` step 36 (§9 item 903, §10 U143): the sphere page's pulsed section is a test.

``docs/guide/sphere_heating_transient.md`` gains "The same heating, pulsed": one
fenced ``python sphere-transient-duty`` block (U135's schedule on U133's study,
Python only), extracted verbatim and exec'd by item 257's helpers (imported from
``test_wf17_first_study_tutorial``, not copied).  Brief:
docs/planning/chunks/WF-17.md, heading "§9 item 903".  Bands imported, none set.

W0 (no solve):
(X)  exactly one tagged block, one ``Study(``, one ``DutyCycle(``, no ``probe``;
(K)  ``Study.solve`` stubbed, ``Study.temperature`` spied (``duty`` recorded from
     ``**kw``): the Study is ``_study(degree=1)``, the cooling is
     ``RobinCooling(K, H, 0.0, source_scale=G)``, both specs are ``EXPECTED_A``,
     call 1 has no duty and call 2 ``DutyCycle(8 * (T_STAR / N_MAIN), 0.5)``;
(E)  the ``sphere-transient-duty-output`` block states ``MAT10_SUB_CELLS`` and
     ``MAT12_A_T_PRINT`` literally;
(D)  docrefs dead 0 and ``MIN_GUIDE_PAGES`` holds.
W1 (``FEM_EM_WF17_S36=1``): (Once), (UT), (χ) binary, (S0) ``_census``, (Δ) the
pulsed θ against Σ_m (−1)^m θ_{H_m} over ten direct runs ≤ ``L_BAND``, (R-A_t)
on the block's ``heating``, (E) ≤ ``E_BAND`` on all twelve runs.
Controls (asserted): (Δ-ctl) dropping m = 9 misses ``L_BAND``; (c-sched) the
pulsed θ misses ``A_BAND`` against the continuous series; (c-page) a removed or
doubled tag raises ``LookupError``; (c-K) ``0.5 -> 1.0`` in the block's
``DutyCycle(`` makes (K)'s mismatch list exactly ``['D.duty.duty']``.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pytest
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver import DutyCycle, RobinCooling, Study, TransientSpec
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import _study
from tests.package.test_wf15_study_sphere_transient import (
    L_BAND,
    MAT12_A_T_PRINT,
    _census,
    _chi_string,
    _on_space,
    _whole_rel,
)
from tests.package.test_wf15_study_temperature import MAT10_SUB_CELLS
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    _study_calls,
    exec_block,
    extract_tagged_block,
)
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
    A_BAND,
    G,
    H,
    K,
    ROBIN_TAG,
    _interior_indicator,
    _interior_rel,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    E_BAND,
    N_MAIN,
    RHO_C,
    T_STAR,
    _exact_theta,
)

TAG = "sphere-transient-duty"
OUTPUT_TAG = "sphere-transient-duty-output"
REPO_ROOT = Path(__file__).resolve().parents[2]
SCRATCH = REPO_ROOT / "logs" / "wf17" / "s36"
S36_ENV = "FEM_EM_WF17_S36"
# U135's schedule (WF-15.md:1260): P = 8 steps of t*/N, D = 0.5.
EXPECTED_DUTY = DutyCycle(8 * (T_STAR / N_MAIN), 0.5)

needs_s36 = pytest.mark.skipif(
    os.environ.get(S36_ENV) != "1", reason=f"W1 solve: set {S36_ENV}=1"
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
    """(K)'s comparator: Study fields, cooling, both specs, and the two duties."""
    ref = _study(degree=1)
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


def test_wf17_s36_extracts_exactly_one_block():
    _, block = _page_and_block()
    calls = _study_calls(block)
    n_duty = block.count("DutyCycle(")
    _say(f"[WF-17s36] (X) Study( calls {len(calls)}; DutyCycle( {n_duty}; "
         f"probe in block: {'probe' in block}")
    assert len(calls) == 1
    assert n_duty == 1
    assert "probe" not in block


def test_wf17_s36_block_configures_u135(monkeypatch):
    _, block = _page_and_block()
    rec = _exec_no_solve(block, monkeypatch)
    mism = _mismatches(rec["study"], rec["calls"])
    duties = [None if d is None else f"DutyCycle({d.period_s:.6g}, {d.duty})"
              for d in (c[2] for c in rec["calls"])]
    _say(f"[WF-17s36] (K) degree={rec['study'].degree} duties={duties} (expected "
         f"[None, DutyCycle({EXPECTED_DUTY.period_s:.6g}, {EXPECTED_DUTY.duty})]) "
         f"mismatches={mism}")
    assert rec["study"].degree == 1
    assert mism == []


def test_wf17_s36_planted_duty_is_caught(monkeypatch):
    _, block = _page_and_block()
    old = "DutyCycle(8 * dt, 0.5)"
    assert block.count(old) == 1
    planted = block.replace(old, "DutyCycle(8 * dt, 1.0)")
    rec = _exec_no_solve(planted, monkeypatch)
    mism = _mismatches(rec["study"], rec["calls"])
    _say(f"[WF-17s36] (c-K) duty 0.5 -> 1.0: mismatches={mism}")
    assert mism == ["D.duty.duty"]


def test_wf17_s36_planted_untagged_and_doubled_pages_raise():
    text, block = _page_and_block()
    untagged = text.replace(f"```python {TAG}\n", "```python\n")
    with pytest.raises(LookupError):
        extract_tagged_block(untagged, TAG)
    doubled = text + f"\n```python {TAG}\n{block}```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    _say("[WF-17s36] (c-page) untagged and doubled pages both raise LookupError")


def test_wf17_s36_expected_output_states_the_records():
    text, _ = _page_and_block()
    out = extract_tagged_block(text, OUTPUT_TAG, lang="text")
    expected = f"sphere cells = {MAT10_SUB_CELLS}\n(A_t) {MAT12_A_T_PRINT}\n"
    _say(f"[WF-17s36] (E) {out!r} vs {expected!r}")
    assert out == expected


def test_wf17_s36_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17s36] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"(MIN_GUIDE_PAGES {MIN_GUIDE_PAGES}) exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ---------------------------------------------------------------- W1 (the solve)


@complex_only
@needs_s36
def test_wf17_s36_block_reproduces_u135(monkeypatch):
    comm = MPI.COMM_WORLD
    _, block = _page_and_block()
    rec: dict = {"calls": [], "thermal": []}
    orig_solve, orig_temp = Study.solve, Study.temperature

    def spy_solve(self, comm=MPI.COMM_WORLD):
        rec["study"] = self
        rec["result"] = orig_solve(self, comm)
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
    _say(f"WF17S36 (K) mismatches={mism}; once {once}; sub cells heating {n_sub} "
         f"pulsed {n_sub_p} (MAT10_SUB_CELLS {MAT10_SUB_CELLS})")

    # (χ) / (S0)
    chis = list(pulsed.source_multipliers)
    chi_binary = all(c in (0.0, 1.0) for c in chis)
    _say(f"WF17S36 (chi) chi string {_chi_string(chis)} (n {len(chis)}; binary {chi_binary}); "
         f"U135 printed its string at 20261108T134930Z_WF-15-s17-W1.log:134")
    _census(pulsed, dt, 2, T_STAR / 2, "WF17S36 pulsed")

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
    _say(f"WF17S36 (Delta) whole rel L2 theta_pulsed - sum(-1)^m theta_H_m = {err_L:.3e} "
         f"(band {L_BAND:.0e}; U135 read 6.608e-16); layouts {lay_H}")
    _say(f"WF17S36 (Delta-ctl) dropping m = 9: {err_L_ctl:.6f} (must miss {L_BAND:.0e}; "
         f"U135 read 0.284785)")

    # (R-A_t) on heating's own submesh; (c-sched) on pulsed's.
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    subA = heating.mesh
    theta_A = _exact_theta(subA, T_STAR, q_cf)
    err_A = _interior_rel(subA, _interior_indicator(subA), heating.temperature - theta_A,
                          theta_A)
    _say(f"WF17S36 (R-A_t) heating interior rel L2 = {err_A:.10e} [{err_A:.4e}] "
         f"(band {A_BAND}; U131 {MAT12_A_T_PRINT})")
    theta_cf = _exact_theta(sub, T_STAR, q_cf)
    err_cs = _interior_rel(sub, _interior_indicator(sub), pulsed.temperature - theta_cf,
                           theta_cf)
    _say(f"WF17S36 (c-sched) theta_pulsed vs continuous series: interior rel L2 = "
         f"{err_cs:.4e} (must miss {A_BAND}; U135 read 4.9989e-01)")

    # Layout check for carrying heating onto pulsed's space (diagnostic only).
    _, lay_A = _on_space(V1, heating.temperature, comm)

    # (E)
    runs = [("heating", heating), ("pulsed", pulsed)] + [
        (f"H{m}", h) for m, h in enumerate(H_runs)]
    es = [(lbl, float(r.max_balance_residual)) for lbl, r in runs]
    _say("WF17S36 (E) max balance " + "  ".join(f"{lbl} {e:.3e}" for lbl, e in es)
         + f"  (band {E_BAND:.0e}; {len(es)} runs)")
    cells_line = f"sphere cells = {MAT10_SUB_CELLS}"
    _say(f"WF17S36 heating/pulsed same layout+geometry {lay_A}; stdout has "
         f"{cells_line!r}: {cells_line in captured}")
    _say(f"WF17S36 elapsed: block (1 EM + 2 transients) {t_block:.1f} s, direct runs (10) "
         f"{t_dir:.1f} s, total {time.perf_counter() - t0:.1f} s")

    assert mism == [], f"(K) {mism}"
    assert once, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS and n_sub_p == MAT10_SUB_CELLS, f"(UT) {n_sub} {n_sub_p}"
    assert chi_binary, f"(chi) {_chi_string(chis)}"
    assert lay_H, "(Delta) layouts differ"
    assert err_L <= L_BAND, f"(Delta) {err_L:.3e} > {L_BAND:.0e}"
    assert err_A <= A_BAND, f"(R-A_t) {err_A:.4e} > {A_BAND}"
    assert f"{err_A:.4e}" == MAT12_A_T_PRINT, f"(R-A_t) print {err_A:.4e}"
    assert all(e <= E_BAND for _, e in es), f"(E) {es}"
    assert err_L_ctl > L_BAND, f"(Delta-ctl) {err_L_ctl:.3e} did not miss {L_BAND:.0e}"
    assert err_cs > A_BAND, f"(c-sched) {err_cs:.4e} did not miss {A_BAND}"
    if comm.rank == 0:
        assert cells_line in captured
    print(f"WF17S36 VERDICT rank {comm.rank} PASS", flush=True)
