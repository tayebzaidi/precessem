"""`WF-1` step 20: ``outputs.probes: {points}`` on item 192's four-port F-small file
drives `POST-8`'s ``point_sar`` (owner-rule default) through ``femem run`` — point
SAR and ``|E|`` per 1 V port drive, the `WF-16` Study's own kept fields.

Brief: docs/planning/chunks/WF-1.md, "§9 item 564" and the review ruling
2026-10-20 03:00.  The solve legs run in the complex build and only with
``FEM_EM_WF1_S20=1`` (then they assert the complex build and the width the
node is for, never skip).  The parse node builds no mesh and runs in either
build.

Anchors: (K) the file's Study = ``f_small_four_port.yaml``'s and the parsed
points by ``float.hex``; (X) refusals naming their field (shape, entries,
outside the phantom, the ball form, capacitor sheets); (R) cells and S tie to
item 192's records; (I-Py) the summary SAR = an independent Python Study's
``point_sar`` at ``L0_BAND``; (U′) the summary is ``point_sar`` on the run's
own fields, bitwise across ranks, JSON = memory, ``|E|`` identity at
``EXACT_IDENTITY_RTOL``; (W) ``-n 1`` / ``-n 4`` vs ``-n 2`` at ``WIDTH_RTOL``.
Control (rule (e), asserted): (c-disp) point 0 displaced 1 mm radially misses
by more than ``WIDTH_RTOL``.  Every band, record and helper is imported, none
restated.
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
    SphereInBoxGeometry,
    StudyConfigError,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.post.sar import point_sar
from fem_em_solver.study_yaml import _load, _study_from_mapping
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf16_study_sar import _port_study
from tests.post.test_point_owner_rule import _bitwise_equal_across_ranks
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_post8_step3a_averaged_row1 import WIDTH_RTOL

S20_ENV = "FEM_EM_WF1_S20"
EM_PROBES = WF1 / "f_small_em_probes.yaml"
FOUR_PORT = WF1 / "f_small_four_port.yaml"
CAPS = WF1 / "f_small_hole_pec_capacitors.yaml"
NPZ = LOGS / "s20" / "em-probes-n2.npz"

# The file's five points (same decimal literals as the YAML).
POINTS = [
    (0.02, 0.0, 0.0),
    (0.0, 0.02, 0.0),
    (-0.02, 0.0, 0.01),
    (0.0, -0.02, -0.01),
    (0.014, 0.014, 0.02),
]
DISPLACED_0 = (0.021, 0.0, 0.0)  # point 0 moved 1 mm radially (c-disp)


def _knob() -> bool:
    raw = os.environ.get(S20_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S20_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _all(comm, flag):
    return bool(comm.allreduce(bool(flag), op=MPI.LAND))


def _field(fn):
    try:
        fn()
    except StudyConfigError as exc:
        return exc.field
    return None


def _hex(points):
    return [[float(x).hex() for x in p] for p in points]


# ------------------------------------------------------------------ no mesh


def test_wf1_s20_parse_and_refusals(monkeypatch, capsys):
    """(K), (X): no mesh, no solve, either build."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (
        SphereInBoxGeometry,
        FHumanBirdcageGeometry,
        FourPortBirdcageGeometry,
        ThermalBallGeometry,
        ThermalBoxGeometry,
    ):
        monkeypatch.setattr(cls, "build", stub)

    diff = _spec_diff(load_study(EM_PROBES), load_study(FOUR_PORT))
    _, kind, outputs, cooling = _load(EM_PROBES)
    got, want = _hex(outputs["probes"]), _hex(POINTS)
    _say(
        capsys,
        f"\n[WF-1 s20] (K) spec diff {diff}; kind {kind!r}; outputs keys "
        f"{sorted(outputs)}; points {outputs['probes']}; hex equal {got == want}",
    )
    assert diff == [], diff
    assert kind == "four_port_birdcage"
    assert cooling is None
    assert set(outputs) == {"probes"}
    assert got == want, (got, want)

    capsys.readouterr()
    rc = cli.main(["check", str(EM_PROBES)])
    lines = capsys.readouterr().out.splitlines()
    rc_all = _rc_all(rc, comm)
    _say(capsys, f"  (K) femem check {EM_PROBES.name}: rc {rc_all}")
    assert rc_all == (0, 0)
    if comm.rank == 0:
        assert [ln for ln in lines if ln.startswith("family: ")] == ["family: four_port_birdcage"]

    s = _mapping(EM_PROBES)
    base = {k: v for k, v in s.items() if k != "outputs"}
    caps = {k: v for k, v in _mapping(CAPS).items() if k != "outputs"}
    good = {"points": [list(p) for p in POINTS]}
    rows = [
        ({**base, "outputs": {"probes": {}}}, "outputs.probes.points"),
        ({**base, "outputs": {"probes": {"points": []}}}, "outputs.probes.points"),
        ({**base, "outputs": {"probes": {"points": [[0.02, 0.0]]}}}, "outputs.probes.points[0]"),
        ({**base, "outputs": {"probes": {**good, "mass": 1.0}}}, "outputs.probes.mass"),
        (
            {**base, "outputs": {"probes": {"points": [[0.02, 0.0, 0.0], [0.0, 0.0, 0.05]]}}},
            "outputs.probes.points[1]",
        ),
        (
            {**base, "outputs": {"probes": {"points": [[0.031, 0.0, 0.0]]}}},
            "outputs.probes.points[0]",
        ),
        (
            {**base, "outputs": {"probes": [{"centre": [0.0, 0.0, 0.0], "radius": 0.002}]}},
            "outputs.probes",
        ),
        ({**caps, "outputs": {"probes": good}}, "outputs.probes"),
        ({**base, "outputs": {"probes": good}}, None),
    ]
    fields = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for (_, want_f), g in zip(rows, fields):
        _say(capsys, f"  (X) want {want_f!r} got {g!r}")
    assert fields == [f for _, f in rows], [(f, g) for (_, f), g in zip(rows, fields) if f != g]
    _say(capsys, f"  build calls {calls['n']}")
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _require(width):
    if not _knob():
        pytest.skip(
            f"{S20_ENV} unset: the four-port probe solve legs are environment-"
            "selected (complex build)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port port Study needs the complex build"
    assert comm.size in width, (comm.size, width)
    return comm


def _arrays(summary):
    order = list(summary["probes"]["values"])
    sar = np.array([summary["probes"]["values"][p]["sar_w_per_kg"] for p in order])
    eab = np.array([summary["probes"]["values"][p]["e_abs_v_per_m"] for p in order])
    return order, sar, eab


def _run(comm, label):
    out = LOGS / label
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    run = run_study_file(EM_PROBES, output_dir=out, comm=comm)
    t_run = time.perf_counter() - t0
    _p0(
        f"\n[WF-1 s20] -n {comm.size}: run_study_file {t_run:.1f} s (summary elapsed "
        f"{run.summary['elapsed_s']:.1f} s); {run.summary['cells']} cells"
    )
    return run


def _legacy(run, comm):
    """``owner=None`` values (printed only; the `POST-8` motivation)."""
    res = run.study.fields_result
    pts = np.asarray(POINTS)
    sigma, rho = float(run.study.material.sigma), float(run.study.density_kg_m3)
    return {
        pid: point_sar(f.e_real, f.e_imag, pts, sigma=sigma, rho=rho, comm=res.comm, owner=None)
        for pid, f in res.fields.items()
    }


def test_wf1_s20_em_probes_n2():
    comm = _require((2,))
    t00 = time.perf_counter()
    run = _run(comm, "em_probes_n2")
    s = run.summary
    pr = s["probes"]
    order, sar, eab = _arrays(s)
    res = run.study.fields_result
    sigma, rho = float(run.study.material.sigma), float(run.study.density_kg_m3)
    pts = np.asarray(POINTS)
    _p0(f"  ports {order}; sigma {sigma!r}; rho {rho!r}; owner {pr['owner']!r}")
    for i, pid in enumerate(order):
        _p0(f"  {pid}: SAR [W/kg per 1 V] {[f'{v:.6e}' for v in sar[i]]}")
        _p0(f"  {pid}: |E| [V/m per 1 V]  {[f'{v:.6e}' for v in eab[i]]}")

    # (R) same fixture
    s_mat = np.array([[complex(re, im) for re, im in row] for row in s["s_matrix"]])
    rel_s = _rel(s_mat, _record(1))
    _p0(
        f"  (R) cells {s['cells']} (record {ANS4_EXACT_STRIP_CELL_RECORD}); "
        f"||S - P1||/||P1|| = {rel_s:.3e} (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    assert s["cells"] == ANS4_EXACT_STRIP_CELL_RECORD
    assert rel_s <= LEG_D0_REPRODUCTION_BAND
    assert [_hex([p])[0] for p in pr["points"]] == _hex(POINTS)
    assert pr["owner"] == "lowest_original_cell"

    # (U′) the summary is point_sar on the run's own fields, same five-point batch
    own = np.array(
        [
            point_sar(
                res.fields[pid].e_real,
                res.fields[pid].e_imag,
                pts,
                sigma=sigma,
                rho=rho,
                comm=res.comm,
            )
            for pid in order
        ]
    )
    exact = bool(np.array_equal(own, sar))
    bitwise = _bitwise_equal_across_ranks(comm, sar)
    js = json.loads(run.paths["summary"].read_text())
    rest = {k: v for k, v in s.items() if k != "elapsed_s"}
    js_rest = {k: v for k, v in js.items() if k != "elapsed_s"}
    same_json = _all(comm, js_rest == rest and (comm.rank != 0 or js == s))
    ident = float(np.max(np.abs(eab**2 * sigma / (2.0 * rho) - sar) / sar))
    _p0(
        f"  (U′) summary == point_sar on run fields {exact}; bitwise across ranks {bitwise}; "
        f"JSON == memory {same_json}; max |E|^2 sigma/(2 rho) vs SAR rel {ident:.3e} "
        f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e})"
    )
    assert exact
    assert bitwise
    assert same_json
    assert ident <= EXACT_IDENTITY_RTOL
    assert run.probes is None
    assert np.all(sar > 0.0)

    # (I-Py) theta_file = theta_Python: an independent Study, its own field solve
    t0 = time.perf_counter()
    twin = _port_study(degree=1)
    tres = twin._solve_port_fields(comm)
    py = np.array(
        [
            point_sar(
                tres.fields[pid].e_real,
                tres.fields[pid].e_imag,
                pts,
                sigma=float(twin.material.sigma),
                rho=float(twin.density_kg_m3),
                comm=tres.comm,
            )
            for pid in order
        ]
    )
    i_py = float(np.max(np.abs(sar - py) / np.abs(py)))
    _p0(
        f"  (I-Py) Python twin {time.perf_counter() - t0:.1f} s, {tres.cells} cells; "
        f"max |SAR_file - SAR_py| / SAR_py = {i_py:.3e} (ASSERTED <= {L0_BAND:.0e})"
    )
    assert i_py <= L0_BAND

    # (c-disp) asserted: point 0 displaced 1 mm radially misses
    disp = np.array(
        [
            point_sar(
                res.fields[pid].e_real,
                res.fields[pid].e_imag,
                np.asarray([DISPLACED_0]),
                sigma=sigma,
                rho=rho,
                comm=res.comm,
            )[0]
            for pid in order
        ]
    )
    c_disp = np.abs(disp - sar[:, 0]) / sar[:, 0]
    _p0(
        f"  (c-disp) 1 mm displaced point 0 per port rel diff {[f'{v:.4e}' for v in c_disp]}; "
        f"min {float(c_disp.min()):.4e} (ASSERTED > {WIDTH_RTOL:.0e}; predicted 1e-2 to 3e-1)"
    )
    assert float(c_disp.min()) > WIDTH_RTOL

    leg = _legacy(run, comm)
    for pid in order:
        _p0(f"  (printed) owner=None {pid}: {[f'{v:.15e}' for v in leg[pid]]}")
        _p0(f"  (printed) owner rule {pid}: " f"{[f'{v:.15e}' for v in sar[order.index(pid)]]}")

    # (W) the -n 2 reference for the width legs: rank 0 writes after an allgather check
    agree = all(
        np.array_equal(a, sar) and np.array_equal(b, eab) for a, b in comm.allgather((sar, eab))
    )
    assert agree
    if comm.rank == 0:
        NPZ.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            NPZ,
            points=np.asarray(POINTS),
            order=np.array(order),
            sar=sar,
            e_abs=eab,
            cells=np.array(s["cells"]),
        )
    comm.Barrier()
    _p0(f"  wrote {NPZ}; node {time.perf_counter() - t00:.1f} s")


def test_wf1_s20_em_probes_width():
    comm = _require((1, 4))
    exists = comm.allreduce(bool(NPZ.exists()), op=MPI.LAND)
    assert (
        exists
    ), f"{NPZ} absent: run test_wf1_s20_em_probes_n2 at -n 2 first (never re-created here)"
    ref = np.load(NPZ)
    t00 = time.perf_counter()
    run = _run(comm, f"em_probes_n{comm.size}")
    s = run.summary
    order, sar, eab = _arrays(s)
    ref_order = [str(x) for x in ref["order"]]
    sp_sar = np.abs(sar / ref["sar"] - 1.0)
    sp_e = np.abs(eab / ref["e_abs"] - 1.0)
    _p0(
        f"  (W) -n {comm.size} vs -n 2: cells {s['cells']} vs {int(ref['cells'])}; port order "
        f"{order} vs {ref_order}"
    )
    for i, pid in enumerate(order):
        _p0(
            f"  (W) {pid}: SAR spread {[f'{v:.3e}' for v in sp_sar[i]]}; "
            f"|E| spread {[f'{v:.3e}' for v in sp_e[i]]}"
        )
    worst = float(max(sp_sar.max(), sp_e.max()))
    _p0(
        f"  (W) max spread {worst:.3e} = {worst / WIDTH_RTOL:.3e} x WIDTH_RTOL "
        f"(ASSERTED <= {WIDTH_RTOL:.0e}); |E| / max|E| per point "
        f"{[f'{v:.3f}' for v in (eab.max(axis=0) / eab.max())]}"
    )
    leg = _legacy(run, comm)
    for pid in order:
        _p0(f"  (printed) owner=None {pid}: {[f'{v:.15e}' for v in leg[pid]]}")
    assert s["cells"] == int(ref["cells"])
    assert order == ref_order
    assert np.array_equal(np.asarray(s["probes"]["points"]), ref["points"])
    assert worst <= WIDTH_RTOL
    _p0(f"  node {time.perf_counter() - t00:.1f} s")
