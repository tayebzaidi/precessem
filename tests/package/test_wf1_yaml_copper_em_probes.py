"""`WF-1` step 22: ``outputs.probes: {points}`` on the copper (Leontovich) F-small hole
file — point SAR and ``|E|`` per 1 V port drive through ``femem run``, on the `ANS-6`
fixture, read by `POST-8`'s ``point_sar`` (owner-rule default).

Brief: docs/planning/chunks/WF-1.md, "§9 item 573" and the review ruling
2026-10-20 10:30 (3).  The solve node runs in the complex build and only with
``FEM_EM_WF1_S22=1`` (then it asserts the complex build and the record width,
never skips).  The parse node builds no mesh and runs in either build.

Anchors: (K) the file's Study = ``f_small_copper.yaml``'s, the points by
``float.hex``, the conductor a ``LeontovichWall`` with the copper file's sigma;
(K1) cells and conductor summary; (R-Cu) the ``cu`` class means at
``LEG_D0_REPRODUCTION_BAND``; (U′) the summary is ``point_sar`` on the run's own
fields, bitwise across ranks, JSON = memory, ``|E|`` identity at
``EXACT_IDENTITY_RTOL``; (I-Py) the summary SAR = an independent Python copper
Study's ``point_sar`` at ``L0_BAND``.  Control (rule (e), asserted): (c-disp)
point 0 displaced 1 mm radially misses by more than ``WIDTH_RTOL``.  Every band,
record and helper is imported, none restated.  No PEC twin (ruling (3)).
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    SphereInBoxGeometry,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.post.sar import point_sar
from fem_em_solver.study_yaml import _load

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_study_copper import _rel_means
from tests.package.test_wf1_study_yaml import LOGS, WF1, _fresh_dir, _rc_all, _spec_diff
from tests.package.test_wf1_yaml_copper import COPPER_FILE, _py_copper
from tests.package.test_wf1_yaml_em_probes import EM_PROBES, POINTS
from tests.post.test_point_owner_rule import _bitwise_equal_across_ranks
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_post8_step3a_averaged_row1 import WIDTH_RTOL

S22_ENV = "FEM_EM_WF1_S22"
CU_PROBES = WF1 / "f_small_copper_em_probes.yaml"
KEY = "64"
DISPLACED_0 = (0.021, 0.0, 0.0)  # point 0 moved 1 mm radially (c-disp)


def _knob() -> bool:
    raw = os.environ.get(S22_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S22_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _all(comm, flag):
    return bool(comm.allreduce(bool(flag), op=MPI.LAND))


def _hex(points):
    return [[float(x).hex() for x in p] for p in points]


# ------------------------------------------------------------------ no mesh


def test_wf1_s22_parse(monkeypatch, capsys):
    """(K): no mesh, no solve, either build."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                ThermalBallGeometry, ThermalBoxGeometry):
        monkeypatch.setattr(cls, "build", stub)

    s_probe, s_cu = load_study(CU_PROBES), load_study(COPPER_FILE)
    diff = _spec_diff(s_probe, s_cu)
    _, kind, outputs, cooling = _load(CU_PROBES)
    _, _, outputs_em, _ = _load(EM_PROBES)
    got, want = _hex(outputs["probes"]), _hex(POINTS)
    c = s_probe.conductor
    _say(capsys, f"\n[WF-1 s22] (K) spec diff {diff}; kind {kind!r}; outputs keys "
         f"{sorted(outputs)}; points {outputs['probes']}; hex equal {got == want}; "
         f"conductor {c!r} sigma hex {c.sigma_s_per_m.hex()} vs copper file "
         f"{s_cu.conductor.sigma_s_per_m.hex()}")
    assert diff == [], diff
    assert kind == "four_port_birdcage"
    assert cooling is None
    assert set(outputs) == {"probes"}
    assert got == want, (got, want)
    assert _hex(outputs_em["probes"]) == got
    assert type(c) is LeontovichWall
    assert c.sigma_s_per_m.hex() == s_cu.conductor.sigma_s_per_m.hex()
    assert s_probe.geometry.as_hole is True

    capsys.readouterr()
    rc = cli.main(["check", str(CU_PROBES)])
    capsys.readouterr()
    rc_all = _rc_all(rc, comm)
    _say(capsys, f"  (K) femem check {CU_PROBES.name}: rc {rc_all}; build calls {calls['n']}")
    assert rc_all == (0, 0)
    assert calls["n"] == 0


# ------------------------------------------------------------------ solve


def test_wf1_s22_copper_em_probes_n2():
    if not _knob():
        pytest.skip(f"{S22_ENV} unset: the copper probe solve is environment-selected "
                    "(complex build)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the copper four-port port Study needs the complex build"
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        comm.size, ANS6_EXACT_STRIP_RECORD_RANK_WIDTH)
    t00 = time.perf_counter()
    out = LOGS / "s22_copper_em_probes_n2"
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    run = run_study_file(CU_PROBES, output_dir=out, comm=comm)
    t_run = time.perf_counter() - t0
    s = run.summary
    pr = s["probes"]
    _p0(f"\n[WF-1 s22] -n {comm.size}: run_study_file {t_run:.1f} s (summary elapsed "
        f"{s['elapsed_s']:.1f} s); {s['cells']} cells; conductor {s.get('conductor')}")

    # (K1) the fixture
    assert s["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert s["conductor"] == {"kind": "leontovich", "sigma_s_per_m": 58000000.0}

    # (R-Cu) the cu class means
    s_mat = np.array([[complex(re, im) for re, im in row] for row in s["s_matrix"]])
    rel = _rel_means(_class_means(s_mat), ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["cu"])
    _p0("  (R-Cu) class means vs cu record, relative: " + "  ".join(
        f"{c} {rel[c]:.3e}" for c in CLASSES) +
        f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}; predicted 4.093e-11 / 8.673e-11 / "
        "1.376e-10)")
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])

    order = list(pr["values"])
    sar = np.array([pr["values"][p]["sar_w_per_kg"] for p in order])
    eab = np.array([pr["values"][p]["e_abs_v_per_m"] for p in order])
    res = run.study.fields_result
    sigma, rho = float(run.study.material.sigma), float(run.study.density_kg_m3)
    pts = np.asarray(POINTS)
    _p0(f"  ports {order}; sigma {sigma!r}; rho {rho!r}; owner {pr['owner']!r}")
    for i, pid in enumerate(order):
        _p0(f"  {pid}: SAR [W/kg per 1 V] {[f'{v:.6e}' for v in sar[i]]}")
        _p0(f"  {pid}: |E| [V/m per 1 V]  {[f'{v:.6e}' for v in eab[i]]}")
    assert [_hex([p])[0] for p in pr["points"]] == _hex(POINTS)
    assert pr["owner"] == "lowest_original_cell"

    # (U′) the summary is point_sar on the run's own fields, same five-point batch
    own = np.array([
        point_sar(res.fields[pid].e_real, res.fields[pid].e_imag, pts, sigma=sigma, rho=rho,
                  comm=res.comm)
        for pid in order
    ])
    exact = bool(np.array_equal(own, sar))
    bitwise = _bitwise_equal_across_ranks(comm, sar)
    js = json.loads(run.paths["summary"].read_text())
    rest = {k: v for k, v in s.items() if k != "elapsed_s"}
    js_rest = {k: v for k, v in js.items() if k != "elapsed_s"}
    same_json = _all(comm, js_rest == rest and (comm.rank != 0 or js == s))
    ident = float(np.max(np.abs(eab ** 2 * sigma / (2.0 * rho) - sar) / sar))
    _p0(f"  (U′) summary == point_sar on run fields {exact}; bitwise across ranks {bitwise}; "
        f"JSON == memory {same_json}; max |E|^2 sigma/(2 rho) vs SAR rel {ident:.3e} "
        f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e})")
    assert exact
    assert bitwise
    assert same_json
    assert ident <= EXACT_IDENTITY_RTOL
    assert np.all(sar > 0.0)

    # (I-Py) theta_file = theta_Python on copper: an independent Study, its own field solve
    t0 = time.perf_counter()
    twin = _py_copper()
    tres = twin._solve_port_fields(comm)
    py = np.array([
        point_sar(tres.fields[pid].e_real, tres.fields[pid].e_imag, pts,
                  sigma=float(twin.material.sigma), rho=float(twin.density_kg_m3),
                  comm=tres.comm)
        for pid in order
    ])
    i_py = float(np.max(np.abs(sar - py) / np.abs(py)))
    _p0(f"  (I-Py) Python copper twin {time.perf_counter() - t0:.1f} s, {tres.cells} cells; "
        f"max |SAR_file - SAR_py| / SAR_py = {i_py:.3e} (ASSERTED <= {L0_BAND:.0e})")
    assert i_py <= L0_BAND

    # (c-disp) asserted: point 0 displaced 1 mm radially misses
    disp = np.array([
        point_sar(res.fields[pid].e_real, res.fields[pid].e_imag, np.asarray([DISPLACED_0]),
                  sigma=sigma, rho=rho, comm=res.comm)[0]
        for pid in order
    ])
    c_disp = np.abs(disp - sar[:, 0]) / sar[:, 0]
    _p0(f"  (c-disp) 1 mm displaced point 0 per port rel diff {[f'{v:.4e}' for v in c_disp]}; "
        f"min {float(c_disp.min()):.4e} (ASSERTED > {WIDTH_RTOL:.0e}; predicted 1e-2 to 3e-1)")
    assert float(c_disp.min()) > WIDTH_RTOL
    _p0(f"  node {time.perf_counter() - t00:.1f} s")
