"""`WF-1` step 28 (§9 item 830, §10 U125): the ``thermal_box`` YAML gains
``thermal.initial`` (θ₀ = g·z, `MAT-21` P47) and ``geometry.wall`` (P48b's
1D-reduced meshed acrylic), and a file reproduces `MAT-21` P63's R1.

Brief: docs/planning/chunks/WF-1.md "§9 item 830" and the review ruling
2026-11-05 03:00 above it.  Real build, ``-n 2``.  The parse node builds no
mesh; the solve node runs only with ``FEM_EM_WF1_S28=1``.

Anchors (bands imported, unmoved): (I) θ at P5 from the R1 file against P63's
direct R1 computed in the same window (P63's own helpers, imported by name, the
same ``READ_PTS`` batch), |Δ|/|θ| ≤ ``IDENTITY_BAND``; the Python ``Study`` the
same; (E) ``max_balance_residual`` ≤ ``E_BAND`` and ``source_power`` equal to
the direct route's; (Cells) 119 808 on both builds; (Ref) the malformed file
refused at ``geometry.wall.thickness_m``, ``femem check`` / ``run`` rc (2, 2).
Controls (asserted): (c-init) the R0 file misses R1 at P5 by > ``IDENTITY_BAND``;
(Drift) that miss equals P63 W1b's printed R1 − R0 at P5, +1.828983e-01 K
(``20261104T060527Z_MAT-21-P63-W1b.log:76``), within 1e-7 (the print's last digit).

Not claimed: the six-wall box, a wall under duty, probes beside a wall, a
gradient on a ball, a cited PMMA ρc (W3).
"""

from __future__ import annotations

import copy
import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    ContainerWall,
    RobinCooling,
    Study,
    StudyConfigError,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    UniformHeatSource,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study import TransverseQuadraticHeatSource
from fem_em_solver.study_yaml import _load, _study_from_mapping
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.probe import ProbeBall
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _cmp,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import RHO_C, T_INF, H, K, Q
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, HALF
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W
from tests.validation.test_mat21_wall_fem import (
    BASE,
    ROBIN,
    _materials,
    _n_wall,
    _robin_facets,
    _walled_box,
)
from tests.validation.test_mat21_wall_pre_equilibration_fem import I_P5, READ_PTS, _read
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR
from tests.validation.test_mat21_wall_transient_fem import RUNG, _rho_c_field

S28_ENV = "FEM_EM_WF1_S28"
R1_FILE = WF1 / "thermal_box_wall_gradient_m4.yaml"
R0_FILE = WF1 / "thermal_box_wall_m4.yaml"
MALFORMED = WF1 / "malformed_thermal_box_wall_without_thickness.yaml"
STEADY = WF1 / "thermal_box_steady_m4.yaml"
DT = T_STAR / N_MAIN
DIV = tuple(b * RUNG for b in BASE)  # (52, 32, 8)
NW = _n_wall(RUNG)  # 2
#: P63 W1b's R1 − R0 at P5 (`20261104T060527Z_MAT-21-P63-W1b.log:76`).
W1B_DRIFT_P5 = 1.828983e-01
DRIFT_BOUND = 1e-7  # the print's last digit
#: P63's cell count (`20261104T060527Z_MAT-21-P63-W1b.log:58`).
CELLS = 119808


def _knob() -> bool:
    raw = os.environ.get(S28_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S28_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _field(fn):
    try:
        fn()
    except StudyConfigError as exc:
        return exc.field
    return None


def _with(m, value, *path):
    m = copy.deepcopy(m)
    d = m
    for p in path[:-1]:
        d = d.setdefault(p, {}) if p not in d else d[p]
    d[path[-1]] = value
    return m


def _without(m, *path):
    m = copy.deepcopy(m)
    d = m
    for p in path[:-1]:
        d = d[p]
    del d[path[-1]]
    return m


def _py_wall():
    return ContainerWall(D_WALL, K_W, RHO_C_W)


def _py_study():
    return Study(
        ThermalBoxGeometry(HALF, DIV, wall=_py_wall(), wall_divisions=NW),
        None,
        UniformHeatSource(Q),
    )


def _raises(exc, fn):
    try:
        fn()
    except exc as e:
        return str(e)
    return None


# ------------------------------------------------------------------ no mesh


def test_wf1_s28_parse_and_refusals(monkeypatch, capsys):
    """W0: literals by float.hex, the spec, the refusals; no mesh, no solve."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    monkeypatch.setattr(ThermalBoxGeometry, "build", stub)
    monkeypatch.setattr(ThermalBallGeometry, "build", stub)

    # premises
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    assert ROBIN == 1

    # (K) literals by float.hex
    m = _mapping(R1_FILE)
    g, t, w = m["geometry"], m["thermal"], m["geometry"]["wall"]
    pairs = {
        "half_widths_m": (
            [float(v).hex() for v in g["half_widths_m"]],
            [float(v).hex() for v in HALF],
        ),
        "divisions": (list(g["divisions"]), list(DIV)),
        "wall.thickness_m": (float(w["thickness_m"]).hex(), D_WALL.hex()),
        "wall.k_w_per_m_k": (float(w["k_w_per_m_k"]).hex(), K_W.hex()),
        "wall.rho_c_j_per_m3_k": (float(w["rho_c_j_per_m3_k"]).hex(), RHO_C_W.hex()),
        "wall.divisions": (w["divisions"], NW),
        "q_w_per_m3": (float(m["excitation"]["q_w_per_m3"]).hex(), Q.hex()),
        "k_w_per_m_k": (float(t["k_w_per_m_k"]).hex(), K.hex()),
        "h_w_per_m2_k": (float(t["h_w_per_m2_k"]).hex(), H_GEL.hex()),
        "t_inf_c": (float(t["t_inf_c"]).hex(), T_INF.hex()),
        "dt_s": (float(t["transient"]["dt_s"]).hex(), DT.hex()),
        "n_steps": (t["transient"]["n_steps"], N_MAIN),
        "rho_c": (float(t["transient"]["rho_c_j_per_m3_k"]).hex(), RHO_C.hex()),
        "theta": (float(t["transient"]["theta"]).hex(), (0.5).hex()),
        "startup": (t["transient"]["startup_half_steps"], 2),
        "gradient": (float(t["initial"]["gradient_k_per_m"]).hex(), G_GRAD.hex()),
    }
    bad = {k: v for k, v in pairs.items() if v[0] != v[1]}
    _say(
        capsys,
        f"\n[WF-1 s28] (K) R1 literals: {len(pairs)} keys, mismatches {bad}; dt "
        f"{float(t['transient']['dt_s']).hex()} == T_STAR/N_MAIN {DT.hex()} "
        f"(repr {DT!r}); half widths (A, B, C) = {(A_BOX, B_BOX, C_BOX)}",
    )
    assert bad == {}, bad
    r0 = _mapping(R0_FILE)
    assert _without(m, "thermal", "initial") == r0  # R0 = R1 without `initial`

    # (K) the Study spec and the fourth element
    for f in (R1_FILE, R0_FILE):
        diff = _spec_diff(load_study(f), _py_study())
        st = load_study(f)
        _say(
            capsys,
            f"  (K) {f.name}: spec diff {diff}; wall {st.geometry.wall}, N_w "
            f"{st.geometry.wall_divisions}",
        )
        assert diff == [], diff
        assert st.geometry.wall_divisions == NW and type(st.geometry.wall_divisions) is int
    want_tr = {
        R1_FILE: TransientSpec(DT, N_MAIN, RHO_C, initial_gradient_k_per_m=G_GRAD),
        R0_FILE: TransientSpec(DT, N_MAIN, RHO_C),
    }
    for f, tr in want_tr.items():
        _, kind, outputs, th = _load(f)
        out = []
        _cmp(th.cooling, RobinCooling(K, H_GEL, T_INF), "cooling", out)
        _cmp(th.transient, tr, "transient", out)
        _say(
            capsys,
            f"  (K) {f.name}: kind {kind!r}, outputs {outputs!r}, duty {th.duty}, "
            f"thermal diffs {out}",
        )
        assert kind == "thermal_box" and outputs == {} and th.duty is None
        assert out == [], out
    # the unwalled box is byte-identical in spec: no wall, no N_w, no gradient
    st = load_study(STEADY)
    assert st.geometry.wall is None and st.geometry.wall_divisions is None

    # (K) femem check on both files
    for f in (R1_FILE, R0_FILE):
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        lines = capsys.readouterr().out.splitlines()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}")
        assert rc_all == (0, 0)
        if comm.rank == 0:
            assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]

    # (Ref) the tracked malformed file
    f_mal = _field(lambda: load_study(MALFORMED))
    d = LOGS / "thermal_box_wall_malformed"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(
        capsys,
        f"  (Ref) {MALFORMED.name}: field {f_mal!r}; run rc {run_all}, check rc "
        f"{check_all}, out exists {exists}",
    )
    assert f_mal == "geometry.wall.thickness_m"
    assert run_all == (2, 2) and check_all == (2, 2) and not exists

    # (Ref) the variant table
    s = _mapping(STEADY)
    rows = [
        # the pinned row (test_wf1_yaml_thermal_box.py:250–251), unedited there
        (
            _with(
                s,
                {"thickness_m": 0.005, "k_w_per_m_k": 0.19, "rho_c_j_per_m3_k": 1.0e6},
                "geometry",
                "wall",
            ),
            "geometry.wall",
        ),
        # no transient: refused at the field before the wall's keys are read
        (_with(s, {"k_w_per_m_k": 0.19}, "geometry", "wall"), "geometry.wall"),
        (_without(m, "geometry", "wall", "divisions"), "geometry.wall.divisions"),
        (_with(m, 0, "geometry", "wall", "divisions"), "geometry.wall.divisions"),
        (_with(m, 2.0, "geometry", "wall", "divisions"), "geometry.wall.divisions"),
        (_with(m, 0.0, "geometry", "wall", "thickness_m"), "geometry.wall.thickness_m"),
        (_with(m, 1.0, "geometry", "wall", "radius"), "geometry.wall.radius"),
        (_with(m, [1, 2], "geometry", "wall"), "geometry.wall"),
        (_with(m, {"period_s": 564.3, "duty": 0.5}, "thermal", "duty"), "thermal.duty"),
        (_with(m, {"directory": "x", "probes": []}, "outputs"), "outputs.probes"),
        (_with(m, {"kind": "uniform"}, "excitation", "source"), "excitation.source"),
        (_with(s, {"gradient_k_per_m": 10.0}, "thermal", "initial"), "thermal.initial"),
        (
            _with(m, "a", "thermal", "initial", "gradient_k_per_m"),
            "thermal.initial.gradient_k_per_m",
        ),
        (
            _with(m, float("inf"), "thermal", "initial", "gradient_k_per_m"),
            "thermal.initial.gradient_k_per_m",
        ),
        (_with(m, 1.0, "thermal", "initial", "offset"), "thermal.initial.offset"),
        (_with(m, {}, "thermal", "initial"), "thermal.initial.gradient_k_per_m"),
    ]
    got = [_field(lambda mm=mm: _study_from_mapping(mm)) for mm, _ in rows]
    for (mm, want_f), gf in zip(rows, got):
        _say(capsys, f"  (Ref) variant -> {gf!r} (want {want_f!r})")
    assert [gf for gf in got] == [wf for _, wf in rows], got

    # (Ref) the API refusals, each before any mesh
    cool = RobinCooling(K, H_GEL, T_INF)
    tr_g = TransientSpec(DT, N_MAIN, RHO_C, initial_gradient_k_per_m=G_GRAD)
    api = {
        "wall without N_w": _raises(
            ValueError, lambda: ThermalBoxGeometry(HALF, DIV, wall=_py_wall())
        ),
        "N_w without wall": _raises(
            ValueError, lambda: ThermalBoxGeometry(HALF, DIV, wall_divisions=2)
        ),
        "N_w = 0": _raises(
            ValueError, lambda: ThermalBoxGeometry(HALF, DIV, wall=_py_wall(), wall_divisions=0)
        ),
        "wall not a ContainerWall": _raises(
            TypeError,
            lambda: ThermalBoxGeometry(HALF, DIV, wall=(0.006, 0.19, 1.0), wall_divisions=2),
        ),
        "gradient not finite": _raises(
            ValueError,
            lambda: TransientSpec(DT, N_MAIN, RHO_C, initial_gradient_k_per_m=float("nan")),
        ),
        "walled steady": _raises(ValueError, lambda: _py_study().temperature(cool)),
        "walled quadratic": _raises(
            ValueError,
            lambda: Study(
                ThermalBoxGeometry(HALF, DIV, wall=_py_wall(), wall_divisions=NW),
                None,
                TransverseQuadraticHeatSource(Q),
            ).temperature(cool, transient=tr_g),
        ),
        "walled probe": _raises(
            ValueError,
            lambda: _py_study().probe([ProbeBall((0.0, 0.0, 0.0), 0.005)], cool, transient=tr_g),
        ),
        "gradient on a ball": _raises(
            ValueError,
            lambda: Study(ThermalBallGeometry(0.05, 0.01), None, UniformHeatSource(Q)).temperature(
                cool, transient=tr_g
            ),
        ),
        "gradient on a quadratic box": _raises(
            ValueError,
            lambda: Study(
                ThermalBoxGeometry(HALF, DIV), None, TransverseQuadraticHeatSource(Q)
            ).temperature(cool, transient=tr_g),
        ),
    }
    for k_, v in api.items():
        _say(capsys, f"  (Ref) API {k_}: {v!r}")
    assert all(v is not None for v in api.values()), api
    _say(capsys, f"  build calls {calls['n']}")
    assert calls["n"] == 0


# ------------------------------------------------------------------ the solve


def test_wf1_s28_wall_gradient_reproduces_p63_r1(capsys):
    """W1: (I), (E), (Cells), (c-init), (Drift) — collected, then asserted."""
    if not _knob():
        pytest.skip(f"{S28_ENV} unset: the walled-box solve is environment-selected")
    comm = MPI.COMM_WORLD
    assert comm.size == 2, comm.size
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    t0 = time.perf_counter()

    # P63's direct R1, its fixture's lines (test_mat21_wall_pre_equilibration_fem.py:202–222)
    tq = time.perf_counter()
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, _inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731
    kw = dict(
        k=kf,
        q=Q * chi,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN,
        theta=0.5,
        startup_half_steps=2,
        dt=DT,
    )
    spec_r1 = dict(rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad)
    direct = solve_transient_conduction(msh, ft, **kw, **spec_r1)
    th_d, ok_d = _read(direct.temperature, READ_PTS)
    cells_d = int(msh.topology.index_map(3).size_global)
    nf_d = comm.allreduce(
        int(np.sum(ft.indices < msh.topology.index_map(2).size_local)), op=MPI.SUM
    )
    _say(
        capsys,
        f"\n[WF-1 s28] direct R1 (P63): cells {cells_d} {tuple(ns)} N_w {nw}; Robin "
        f"facets {nf_d}; (E) {direct.max_balance_residual:.3e}; source_power "
        f"{float(direct.source_power):.12e} W; theta(P5) {th_d[I_P5]:+.12e} K; "
        f"{time.perf_counter() - tq:.1f} s",
    )

    runs = {}
    for key, f in (("R1", R1_FILE), ("R0", R0_FILE)):
        tq = time.perf_counter()
        out = LOGS / f"thermal_box_wall_{key}"
        _fresh_dir(out, comm)
        run = run_study_file(f, output_dir=out, comm=comm)
        th, ok = _read(run.thermal.temperature, READ_PTS)
        on_disk = None
        if comm.rank == 0:
            on_disk = json.loads((out / "summary.json").read_text())
        on_disk = comm.bcast(on_disk, root=0)
        runs[key] = dict(run=run, th=th, ok=ok, disk=on_disk, s=time.perf_counter() - tq)
        tm = run.summary["temperature"]
        _say(
            capsys,
            f"  file {f.name}: cells {run.summary['cells']}; (E) "
            f"{tm['max_balance_residual']:.3e}; source_power {tm['source_power_w']:.12e} W; "
            f"theta(P5) {th[I_P5]:+.12e} K; wall {tm.get('wall')}; gradient "
            f"{tm.get('initial_gradient_k_per_m')}; {runs[key]['s']:.1f} s",
        )

    tq = time.perf_counter()
    py = _py_study().temperature(
        RobinCooling(K, H_GEL, T_INF),
        transient=TransientSpec(DT, N_MAIN, RHO_C, initial_gradient_k_per_m=G_GRAD),
    )
    th_py, ok_py = _read(py.temperature, READ_PTS)
    _say(
        capsys,
        f"  Python Study R1: theta(P5) {th_py[I_P5]:+.12e} K; (E) "
        f"{py.max_balance_residual:.3e}; {time.perf_counter() - tq:.1f} s",
    )

    r1, r0 = runs["R1"], runs["R0"]
    ref = th_d[I_P5]
    rel_i = abs(r1["th"][I_P5] - ref) / abs(ref)
    rel_py = abs(th_py[I_P5] - ref) / abs(ref)
    rel_r0 = abs(r0["th"][I_P5] - ref) / abs(ref)
    drift = ref - r0["th"][I_P5]
    src_d = float(direct.source_power)
    src_rel = {
        k: abs(v["run"].summary["temperature"]["source_power_w"] - src_d) / src_d
        for k, v in runs.items()
    }
    batch_rel = float(np.max(np.abs(r1["th"] - th_d) / np.abs(th_d).clip(min=1e-300)))
    _say(
        capsys,
        f"  (I) file R1 vs direct at P5 rel {rel_i:.3e}; Python vs direct {rel_py:.3e} "
        f"(IDENTITY_BAND {IDENTITY_BAND:.0e}); whole READ_PTS batch max rel {batch_rel:.3e} "
        f"(printed)",
    )
    _say(
        capsys,
        f"  (E) source_power rel vs direct {src_rel}; (c-init) R0 file vs R1 at P5 rel "
        f"{rel_r0:.3e} (asserted > {IDENTITY_BAND:.0e}); (Drift) R1 - R0 at P5 {drift:+.9e} K "
        f"vs W1b {W1B_DRIFT_P5:+.6e}, |delta| {abs(drift - W1B_DRIFT_P5):.3e} (bound "
        f"{DRIFT_BOUND:g})",
    )
    if rel_i > IDENTITY_BAND:  # the brief's per-ingredient print on an (I) red
        g1 = r1["run"].thermal
        _say(
            capsys,
            f"  INGREDIENTS file: cells {g1.mesh.topology.index_map(3).size_global}, "
            f"Robin facets {comm.allreduce(int(np.sum(g1.facet_tags.indices < g1.mesh.topology.index_map(2).size_local)), op=MPI.SUM)}, "
            f"sum k {comm.allreduce(float(np.sum(g1.k.x.array[:g1.mesh.topology.index_map(3).size_local])), op=MPI.SUM):.12e} "
            f"vs direct {comm.allreduce(float(np.sum(kf.x.array[:msh.topology.index_map(3).size_local])), op=MPI.SUM):.12e}, "
            f"sum rho c {comm.allreduce(float(np.sum(g1.rho_c.x.array[:g1.mesh.topology.index_map(3).size_local])), op=MPI.SUM):.12e} "
            f"vs direct {comm.allreduce(float(np.sum(rcf.x.array[:msh.topology.index_map(3).size_local])), op=MPI.SUM):.12e}, "
            f"theta0(P5) {G_GRAD * READ_PTS[I_P5][2]:+.6e}",
        )

    disk_diffs = {}
    for k_, v in runs.items():
        mine = json.loads(json.dumps(v["run"].summary))
        disk_diffs[k_] = sorted(
            kk
            for kk in set(mine) | set(v["disk"])
            if kk != "elapsed_s" and mine.get(kk) != v["disk"].get(kk)
        )
    all_diffs = comm.gather(disk_diffs, root=0)
    _say(
        capsys,
        f"  summary vs summary.json differing keys (elapsed_s excluded) per rank " f"{all_diffs}",
    )
    gates = {
        "valid_points": all(ok_d) and all(r1["ok"]) and all(r0["ok"]) and all(ok_py),
        "I_file": rel_i <= IDENTITY_BAND,
        "I_python": rel_py <= IDENTITY_BAND,
        "E_direct": float(direct.max_balance_residual) <= E_BAND,
        "E_python": float(py.max_balance_residual) <= E_BAND,
        "Cells": cells_d == CELLS
        and all(v["run"].summary["cells"] == CELLS for v in runs.values()),
        # rank 0 writes the file; ``elapsed_s`` is each rank's own (first W1,
        # `20261105T210601Z_WF-1-s28-W1.log`: rank 1 red on this gate only; item 390b's
        # record, ranks >= 1 differ in ``elapsed_s`` only), so it is left out of the identity.
        "summary_on_disk": all(diffs == [] for diffs in disk_diffs.values()),
        "summary_keys": "initial_gradient_k_per_m" in r1["run"].summary["temperature"]
        and "initial_gradient_k_per_m" not in r0["run"].summary["temperature"],
        "c_init_misses": rel_r0 > IDENTITY_BAND,
        "Drift": abs(drift - W1B_DRIFT_P5) <= DRIFT_BOUND,
    }
    for k_, v in runs.items():
        tm = v["run"].summary["temperature"]
        gates[f"E_{k_}"] = tm["max_balance_residual"] <= E_BAND
        gates[f"source_power_{k_}"] = src_rel[k_] <= IDENTITY_BAND
    _say(capsys, f"  gates {gates}; elapsed {time.perf_counter() - t0:.1f} s")
    assert all(gates.values()), gates
