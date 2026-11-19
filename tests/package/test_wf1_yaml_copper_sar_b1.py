"""`WF-1` step 8: copper 10 g SAR and B1+ through the YAML front end, and
solid-route refusal parity, answering the untested YAML conductor plus outputs
combination named in `WF-16` item 388's scope.

Brief: docs/planning/chunks/WF-1.md, "§9 item 394" (chain U44, provisional).
``run_study_file`` on ``f_small_copper_sar_b1.yaml`` (the F-small exact-
rectangle copper hole Study, 64 MHz, degree 1, ``-n 2``) with ``outputs.sar`` /
``outputs.b1_plus`` at the origin.  Asserted:

(K)    the file parses bitwise to the Python copper Study and to
       ``f_small_copper.yaml``'s Study; the outputs parse as written;
(X=)   a solid-route ``leontovich`` file refuses at ``conductor`` with a
       ``ValueError`` naming ``as_hole`` (cli rc 2, nothing written, no mesh),
       as ``Study.sar`` / ``b1_plus_at`` / ``solve()`` refuse on the Python side
       (parity of outcome, not of call site); its control flips that one bit;
(K1)   the YAML Study's field route carries the conductor (cells, one term,
       facet tags, summary block);
(S=Y)  its field-route S equals its ``solve()`` S at ``LEG_D0_REPRODUCTION_BAND``;
       (R) the field route reproduces the ``cu`` class means at the same band;
(I-SAR)/(I-B) the summary's SAR / ``|B1+|`` equal the Python copper Study's
       field route per drive at ``REUSE_REPRODUCTION_RTOL``; (iii)/(B2') equal
       the API's values on the same Study at ``EXACT_IDENTITY_RTOL``;
(P)    the kept field's copper power identity at ``DISCRETE_IDENTITY_RTOL``.
Negative controls (rule (e), asserted, backed): (N1) the PEC-hole Study's SAR
and ``|B1+|`` at P1 / origin miss the YAML Cu values by > ``REUSE_REPRODUCTION_RTOL``;
(N2) the P_surf-blind residual is >= ``BLIND_FACTOR`` x ``DISCRETE_IDENTITY_RTOL``.
Item 388's printed SAR / ``|B1+|`` records are printed beside, never asserted.
Every band is imported; none is set here.
"""

from __future__ import annotations

import copy
import dataclasses
import json
import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
import fem_em_solver.study_yaml as study_yaml
from fem_em_solver import (
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    Study,
    StudyConfigError,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study_yaml import _study_from_mapping
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_copper import COPPER_FILE, _py_copper
from tests.package.test_wf16_study_conductor_fields import CENTRE, _kept_identity, _phantom
from tests.package.test_wf16_study_copper import _geometry, _rel_means
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR

KEY = "64"
CU_SAR_FILE = WF1 / "f_small_copper_sar_b1.yaml"
MALFORMED = WF1 / "malformed_conductor_solid_with_outputs.yaml"
SAR_KEY = "averaged_sar_w_per_kg"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def test_wf1_s8_parse_and_refusal_parity(monkeypatch):
    """(K), (X=) and its control: no mesh on either side."""
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=None):
        calls["n"] += 1
        return {"mesh": None}

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", stub)

    # (K) bitwise up to the solve.
    s_yaml = load_study(CU_SAR_FILE)
    d_py = _spec_diff(s_yaml, _py_copper())
    d_cu = _spec_diff(s_yaml, load_study(COPPER_FILE))
    _, _, outputs, _ = study_yaml._load(CU_SAR_FILE)
    _p0(
        f"\n[WF-1 s8] (K) diff vs Python Cu {d_py}; diff vs f_small_copper.yaml {d_cu}; "
        f"outputs {outputs}; mass_g hex {outputs['sar']['mass_g'].hex()}"
    )
    assert d_py == [], d_py
    assert d_cu == [], d_cu
    assert outputs["sar"] == {"mass_g": 10.0, "centres": [(0.0, 0.0, 0.0)]}, outputs
    assert outputs["sar"]["mass_g"].hex() == (10.0).hex()
    assert outputs["b1_plus"]["points"] == [(0.0, 0.0, 0.0)], outputs
    assert set(outputs) == {"sar", "b1_plus"}

    # (X=) YAML side.
    with pytest.raises(StudyConfigError) as ei:
        load_study(MALFORMED)
    exc = ei.value
    _p0(f"  (X=) YAML: type {type(exc).__name__}, field {exc.field!r}; {exc}")
    assert exc.field == "conductor"
    assert isinstance(exc, ValueError)
    assert "as_hole" in str(exc)
    assert calls["n"] == 0
    out = LOGS / "malformed_conductor_solid_with_outputs"
    _fresh_dir(out, comm)
    rc_run = cli.main(["run", str(MALFORMED), "--output-dir", str(out)])
    rc_check = cli.main(["check", str(MALFORMED)])
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
    n_build = comm.allreduce(calls["n"], op=MPI.MAX)
    _p0(
        f"  (X=) cli run rc {run_all}, check rc {check_all}, out exists {exists}, "
        f"build calls {n_build}"
    )
    assert run_all == (2, 2)
    assert check_all == (2, 2)
    assert not exists
    assert n_build == 0

    # (X=)'s control: that one bit alone; and the refusal is the conductor's.
    m = copy.deepcopy(_mapping(MALFORMED))
    assert m["geometry"]["as_hole"] is False
    m["geometry"]["as_hole"] = True
    d_ctl = _spec_diff(_study_from_mapping(m, "<as_hole true>"), s_yaml)
    m2 = copy.deepcopy(_mapping(MALFORMED))
    del m2["outputs"]
    with pytest.raises(StudyConfigError) as ei2:
        _study_from_mapping(m2, "<no outputs>")
    _p0(
        f"  (X=) control: as_hole true diff {d_ctl}; outputs deleted -> field "
        f"{ei2.value.field!r}"
    )
    assert d_ctl == [], d_ctl
    assert ei2.value.field == "conductor"
    assert "as_hole" in str(ei2.value)

    # (X=) Python side: the same combination refuses at its first solve-bearing call.
    solid = Study(
        dataclasses.replace(_geometry(), as_hole=False),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        conductor=LeontovichWall(study_mod.COPPER_SIGMA),
    )
    with pytest.raises(ValueError, match="as_hole") as e1:
        solid.sar(mass_g=10.0, centres=[CENTRE])
    assert getattr(solid, "fields_result", None) is None
    with pytest.raises(ValueError, match="as_hole"):
        solid.b1_plus_at([CENTRE])
    assert getattr(solid, "fields_result", None) is None
    with pytest.raises(ValueError, match="as_hole"):
        solid.solve()
    assert getattr(solid, "fields_result", None) is None
    assert getattr(solid, "result", None) is None
    # Any build call here is the stub's (no mesh); the count is printed, not asserted.
    _p0(
        f"  (X=) Python: sar / b1_plus_at / solve() refuse ({e1.value}); stub build calls "
        f"{calls['n']} (each returned no mesh)"
    )


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"(R) is asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    assert th14.FREQUENCIES_HZ[KEY] == FREQUENCY_64_HZ
    out = {}
    # (1) the YAML run: solve() then the field route (two builds, ruling 2026-10-08 (4)).
    d = LOGS / CU_SAR_FILE.stem
    _fresh_dir(d, comm)
    t0 = time.perf_counter()
    run = run_study_file(CU_SAR_FILE, output_dir=d, comm=comm)
    out["t_run"] = time.perf_counter() - t0
    with open(d / "summary.json") as fh:
        summary = json.load(fh)
    fr = run.study.fields_result
    _p0(
        f"\n[WF-1 s8] run_study_file {out['t_run']:.1f} s (summary elapsed_s "
        f"{summary['elapsed_s']:.1f}); solve() cells {run.result.cells}; field route cells "
        f"{fr.cells} (mesh {fr.mesh_time:.1f} s, sweep {fr.sweep_time:.1f} s), -n {comm.size}"
    )
    # (2) the same-Study API re-read (no third build).
    t0 = time.perf_counter()
    api_sar = run.study.sar(mass_g=10.0, centres=[CENTRE], comm=comm)
    api_b1 = run.study.b1_plus_at([CENTRE], comm=comm)
    out["api_same"] = run.study.fields_result is fr
    _p0(
        f"[WF-1 s8] same-Study API re-read {time.perf_counter() - t0:.1f} s; "
        f"fields_result reused {out['api_same']}"
    )
    # (3) Python Cu Study A, field route only.
    a = _py_copper()
    t0 = time.perf_counter()
    py_sar = a.sar(mass_g=10.0, centres=[CENTRE], comm=comm)
    py_b1 = a.b1_plus_at([CENTRE], comm=comm)
    fa = a.fields_result
    _p0(
        f"[WF-1 s8] Python Cu Study A field route + SAR/B1 {time.perf_counter() - t0:.1f} s "
        f"(mesh {fa.mesh_time:.1f}, sweep {fa.sweep_time:.1f})"
    )
    # (4) Python PEC-hole Study B, P1 only.
    b = _py_copper("pec")
    t0 = time.perf_counter()
    pec_sar = b.sar(mass_g=10.0, centres=[CENTRE], drives=["P1"], comm=comm)
    pec_b1 = b.b1_plus_at([CENTRE], drives=["P1"], comm=comm)
    fb = b.fields_result
    _p0(
        f"[WF-1 s8] Python PEC Study B field route + SAR/B1 {time.perf_counter() - t0:.1f} s "
        f"(mesh {fb.mesh_time:.1f}, sweep {fb.sweep_time:.1f})"
    )
    # (5) the kept-field identity on the YAML Study's field.
    t0 = time.perf_counter()
    power = _kept_identity(fr, run.study.omega, comm)
    _p0(
        f"[WF-1 s8] kept-field identity {time.perf_counter() - t0:.1f} s: "
        f"P_src {power['p_src']:.9e} W, P_sheets {power['p_sheets']:.9e} W, "
        f"P_phantom {power['p_phantom']:.9e} W, P_non_phantom {power['p_non_phantom']:.9e} W, "
        f"P_surf {power['p_surf']:.9e} W, residual {power['residual']:.3e} "
        f"(band {DISCRETE_IDENTITY_RTOL:.0e}), blind residual {power['blind_residual']:.4e}"
    )
    _p0(
        "[WF-1 s8] prior (Python kept field, 20261012T033238Z_WF-16-s11-W1.log:2692): residual "
        "8.505e-13, P_src 4.806472101e-03 W, P_surf 1.516873149e-06 W, "
        "P_phantom 1.870050049e-06 W"
    )
    ysar = summary["sar"]["values"]["P1"][0]
    yb1 = summary["b1_plus"]["values"]["P1"][0]
    _p0(
        f"[WF-1 s8] printed only: YAML SAR10g(P1, origin) {ysar[SAR_KEY]:.9e} W/kg beside item "
        f"388's Cu record 2.408073397e-06; YAML |B1+|(P1, origin) {yb1:.9e} T beside "
        f"1.833178311e-08 (per 1 V)"
    )
    out.update(
        run=run,
        summary=summary,
        fr=fr,
        api_sar=api_sar,
        api_b1=api_b1,
        py_sar=py_sar,
        py_b1=py_b1,
        pec_sar=pec_sar,
        pec_b1=pec_b1,
        power=power,
    )
    return out


@complex_only
def test_wf1_s8_yaml_field_route_applies_conductor(runs):
    """(K1), (S=Y), (R)."""
    fr, run, summary = runs["fr"], runs["run"], runs["summary"]
    assert fr.cells == summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert fr.conductor_terms is not None and len(fr.conductor_terms) == 1
    assert fr.conductor_facet_tags is not None
    assert summary["conductor"] == {"kind": "leontovich", "sigma_s_per_m": study_mod.COPPER_SIGMA}
    s_f = np.asarray(fr.s_matrix)
    s_s = np.asarray(run.result.s_matrix)
    rel = float(np.linalg.norm(s_f - s_s) / np.linalg.norm(s_s))
    means = _rel_means(_class_means(fr.s_matrix), ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    _p0(
        f"\n[WF-1 s8] (S=Y) ||S_field - S_solve||_F / ||S_solve||_F {rel:.3e} (band "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}; Python prior 7.820e-16); (R) field route vs cu "
        "record: " + "  ".join(f"{c} {means[c]:.3e}" for c in CLASSES)
    )
    assert rel <= LEG_D0_REPRODUCTION_BAND, rel
    for c in CLASSES:
        assert means[c] <= LEG_D0_REPRODUCTION_BAND, (c, means[c])


@complex_only
def test_wf1_s8_yaml_sar_reproduces_python(runs):
    """(I-SAR) and (iii)."""
    summary, run = runs["summary"], runs["run"]
    ys = summary["sar"]["values"]
    py, api = runs["py_sar"], runs["api_sar"]
    assert list(ys) == list(py) == ["P1", "P2", "P3", "P4"], (list(ys), list(py))
    worst, worst_iii = 0.0, 0.0
    for pid in ys:
        assert len(ys[pid]) == len(py[pid]) == 1
        y, p, a = ys[pid][0], py[pid][0], api[pid][0]
        keys = [k for k, v in p.items() if float(v) != 0.0]
        assert SAR_KEY in keys
        rels = {k: abs(float(y[k]) - float(p[k])) / abs(float(p[k])) for k in keys}
        rel_iii = {k: abs(float(y[k]) - float(a[k])) / abs(float(a[k])) for k in keys}
        _p0(
            f"\n[WF-1 s8] (I-SAR) {pid}: SAR10g YAML {float(y[SAR_KEY]):.9e} Python "
            f"{float(p[SAR_KEY]):.9e}; max rel over {len(keys)} keys {max(rels.values()):.3e} "
            f"(band {REUSE_REPRODUCTION_RTOL:.0e}); (iii) vs same-Study API "
            f"{max(rel_iii.values()):.3e} (band {EXACT_IDENTITY_RTOL:.0e})"
        )
        worst = max(worst, max(rels.values()))
        worst_iii = max(worst_iii, max(rel_iii.values()))
        for k in keys:
            assert rels[k] <= REUSE_REPRODUCTION_RTOL, (pid, k, rels[k])
            assert rel_iii[k] <= EXACT_IDENTITY_RTOL, (pid, k, rel_iii[k])
    _p0(f"[WF-1 s8] (I-SAR) max over four drives {worst:.3e}; (iii) max {worst_iii:.3e}")
    assert summary["sar"] == run.summary["sar"]
    assert summary["b1_plus"] == run.summary["b1_plus"]


@complex_only
def test_wf1_s8_yaml_b1_reproduces_python(runs):
    """(I-B) and (B2')."""
    yb = runs["summary"]["b1_plus"]
    py, api = runs["py_b1"], runs["api_b1"]
    assert list(yb["values"]) == list(py) == ["P1", "P2", "P3", "P4"]
    worst = 0.0
    for pid in yb["values"]:
        y = float(yb["values"][pid][0])
        p = float(py[pid][0][0])
        rel = abs(y - p) / abs(p)
        worst = max(worst, rel)
        _p0(
            f"\n[WF-1 s8] (I-B) {pid}: |B1+| YAML {y:.9e} Python {p:.9e} rel {rel:.3e} "
            f"(band {REUSE_REPRODUCTION_RTOL:.0e})"
        )
        assert yb["valid"][pid] == [True], (pid, yb["valid"][pid])
        assert bool(py[pid][1][0]), pid
        assert rel <= REUSE_REPRODUCTION_RTOL, (pid, rel)
    y_all = np.array([float(yb["values"][pid][0]) for pid in yb["values"]])
    a_all = np.array([float(api[pid][0][0]) for pid in yb["values"]])
    b2 = float(np.max(np.abs(y_all - a_all)) / np.max(np.abs(a_all)))
    _p0(
        f"[WF-1 s8] (I-B) max over four drives {worst:.3e}; (B2') vs same-Study API "
        f"{b2:.3e} (band {EXACT_IDENTITY_RTOL:.0e})"
    )
    assert b2 <= EXACT_IDENTITY_RTOL, b2


@complex_only
def test_wf1_s8_kept_field_power_identity(runs):
    """(P) on the YAML Study's kept field."""
    p = runs["power"]
    assert runs["api_same"] is True
    assert p["p_surf"] > 0.0, p
    assert p["residual"] <= DISCRETE_IDENTITY_RTOL, p


@complex_only
def test_wf1_s8_negative_pec_column(runs):
    """(N1) asserted, backed: the PEC-hole Study misses the YAML Cu values."""
    ysar = float(runs["summary"]["sar"]["values"]["P1"][0][SAR_KEY])
    yb1 = float(runs["summary"]["b1_plus"]["values"]["P1"][0])
    psar = float(runs["pec_sar"]["P1"][0][SAR_KEY])
    pb1 = float(runs["pec_b1"]["P1"][0][0])
    r_sar = abs(ysar - psar) / abs(psar)
    r_b1 = abs(yb1 - pb1) / abs(pb1)
    _p0(
        f"\n[WF-1 s8] (N1) PEC SAR10g(P1, origin) {psar:.9e} vs YAML Cu {ysar:.9e}: rel "
        f"{r_sar:.3e} = {r_sar / REUSE_REPRODUCTION_RTOL:.1f}x band; PEC |B1+| {pb1:.9e} vs "
        f"{yb1:.9e}: rel {r_b1:.3e} = {r_b1 / REUSE_REPRODUCTION_RTOL:.1f}x band "
        "(backing 4.151e-04 / 3.578e-04)"
    )
    assert r_sar > REUSE_REPRODUCTION_RTOL, r_sar
    assert r_b1 > REUSE_REPRODUCTION_RTOL, r_b1


@complex_only
def test_wf1_s8_negative_surface_blind(runs):
    """(N2) asserted, backed: the P_surf-blind residual clears the 10x bar."""
    p = runs["power"]
    bar = BLIND_FACTOR * DISCRETE_IDENTITY_RTOL
    _p0(
        f"\n[WF-1 s8] (N2) blind residual {p['blind_residual']:.4e} = "
        f"{p['blind_residual'] / DISCRETE_IDENTITY_RTOL:.1f}x DISCRETE_IDENTITY_RTOL "
        f"(bar {bar:.0e}; backing 3.1559e-04)"
    )
    assert p["blind_residual"] >= bar, p
