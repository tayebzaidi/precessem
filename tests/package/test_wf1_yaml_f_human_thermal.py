"""`WF-1` step 10a (§9 item 431a): `MAT-15` / `MAT-20`'s F-human heating, schema only.

A ``thermal:`` section on an ``f_human_birdcage`` file now parses into a
transient Robin-cooled temperature of the tuned copper F-human's phantom:
the cooling, a ``transient`` grid, 10 g ``balls`` and an optional
``window_s`` / ``sense`` (brief: docs/planning/chunks/WF-1.md, "§9 item
431", part 431a).  ``run_study_file``'s F-human branch calls
``Study.temperature`` on the kept field and writes the per-step series into
the summary.  Nothing here builds a mesh or solves.  Asserted:

(K)     ``load_study`` on the tracked thermal file equals the Python ``Study``
        and the 390a file field by field; the parsed ``RobinCooling``,
        ``TransientSpec``, ball centre, mass and window equal the P17 / step-1
        / MAT-12 / MAT-20 constants and the series JSON by ``float.hex``;
        ``cli check`` returns 0 and prints the family;
(X)     one-key mutations raise ``StudyConfigError`` naming the field; run /
        check on a refused file return 2 on every rank, nothing written;
(Run-0) the branch on counting fakes: one tuned solve, one temperature call,
        the summary equals the fake's values, the JSON equals the summary;
        without ``thermal:`` no temperature call and no ``temperature`` key.
Negative controls (asserted, deterministic): (c1) a dropped ``source_scale``
reads ``RobinCooling``'s default 1.0; (c2) a one-ULP ``source_scale`` differs
by exactly one ULP.

The series JSON is gitignored: read on rank 0 and broadcast; absent ⇒ FAIL.
"""

from __future__ import annotations

import copy
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import RobinCooling, StudyConfigError, TwoPortQuadratureDrive, cli, load_study
from fem_em_solver import study as study_mod
from fem_em_solver.study import TransientSpec
from fem_em_solver.study_yaml import _load, _study_from_mapping, run_study_file
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _p0,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_f_human import FH_FILE
from tests.package.test_wf16_study_tuned_f_human import _study
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    T_STAR,
)
from tests.validation.test_mat15_f_human_tuned_heating import SERIES_JSON
from tests.validation.test_mat20_running_average import W_B2

FHT_FILE = WF1 / "f_human_tuned_copper_thermal.yaml"
REPO = Path(__file__).resolve().parents[2]


def _yf(x) -> str:
    r = repr(float(x))
    if "." not in r and "inf" not in r and "nan" not in r:
        r = "%.17e" % float(x)
    return r


def _series_json():
    """The P17 series JSON, read on rank 0 and broadcast (None if absent)."""
    comm = MPI.COMM_WORLD
    data = None
    if comm.rank == 0:
        path = REPO / SERIES_JSON
        if path.is_file():
            data = json.loads(path.read_text())
    return comm.bcast(data, root=0)


def _literal_lines() -> str:
    """The YAML text of the ``thermal:`` literals (printed, never typed)."""
    js = _series_json()
    assert js is not None, f"{SERIES_JSON} is absent (gitignored; item 380 writes it)"
    c = js["c_star_m"]
    lines = [
        f"# source: {js['source']}",
        "thermal:",
        f"  k_w_per_m_k: {_yf(K_GEL)}",
        f"  h_w_per_m2_k: {_yf(H_GEL)}",
        f"  t_inf_c: {_yf(T_INF)}",
        f"  source_scale: {_yf(js['G_per_w'])}",
        "  transient:",
        f"    dt_s: {_yf(T_STAR / N_MAIN)}",
        f"    n_steps: {int(N_MAIN)}",
        f"    rho_c_j_per_m3_k: {_yf(RHO_C)}",
        "    theta: 0.5",
        "    startup_half_steps: 2",
        "  balls:",
        "    mass_g: 10.0",
        "    centres:",
        f"      - [{', '.join(_yf(v) for v in c)}]",
        f"  window_s: {_yf(W_B2)}",
    ]
    return "\n".join(lines)


def _hex(x):
    return float(x).hex()


def _raises(mapping, source="<mutation>"):
    with pytest.raises(StudyConfigError) as ei:
        _study_from_mapping(mapping, source)
    return ei.value


def test_wf1_fht_parse_and_refusals():
    """(K), (X), (c1), (c2) — no mesh."""
    import yaml

    comm = MPI.COMM_WORLD
    js = _series_json()
    _p0(f"\n[WF-1 s10a] series JSON {SERIES_JSON} present: {js is not None}")
    assert js is not None, "the series JSON is absent: FAIL, never skip"
    _p0("[WF-1 s10a] expected literals:\n" + _literal_lines())

    # ---- (K) bitwise up to the solve --------------------------------------
    s_yaml = load_study(FHT_FILE)
    d_py, d_fh = _spec_diff(s_yaml, _study()), _spec_diff(s_yaml, load_study(FH_FILE))
    _p0(f"  (K) diff vs Python Study {d_py}; vs 390a file {d_fh}")
    assert d_py == [] and d_fh == []
    _, kind, outputs, fh = _load(FHT_FILE)
    assert kind == "f_human_birdcage" and outputs == {}
    exp_cool = RobinCooling(K_GEL, H_GEL, T_INF, source_scale=js["G_per_w"])
    exp_tr = TransientSpec(
        dt_s=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        rho_c_j_per_m3_k=RHO_C,
        theta=0.5,
        startup_half_steps=2,
    )
    cool, tr = fh.cooling, fh.transient
    cool_hex = {
        k: (_hex(getattr(cool, k)), _hex(getattr(exp_cool, k)))
        for k in ("k_w_per_m_k", "h_w_per_m2_k", "t_inf_c", "source_scale")
    }
    tr_hex = {
        k: (_hex(getattr(tr, k)), _hex(getattr(exp_tr, k)))
        for k in ("dt_s", "rho_c_j_per_m3_k", "theta")
    }
    for k, (g, e) in {**cool_hex, **tr_hex}.items():
        _p0(f"  (K) {k}: parsed {g} expected {e}")
    assert type(cool) is RobinCooling and type(tr) is TransientSpec
    assert all(g == e for g, e in cool_hex.values())
    assert cool.degree == exp_cool.degree
    assert all(g == e for g, e in tr_hex.values())
    assert (tr.n_steps, tr.startup_half_steps, tr.t_init_c) == (N_MAIN, 2, None)
    got_c = fh.balls["centres"]
    _p0(
        f"  (K) ball centres {[[_hex(v) for v in c] for c in got_c]} vs JSON "
        f"{[_hex(v) for v in js['c_star_m']]}; mass_g {_hex(fh.balls['mass_g'])}; window_s "
        f"{_hex(fh.window_s)} vs W_B2 {_hex(W_B2)}; sense {fh.sense!r}"
    )
    assert len(got_c) == 1
    assert [_hex(v) for v in got_c[0]] == [_hex(v) for v in js["c_star_m"]]
    assert _hex(fh.balls["mass_g"]) == _hex(10.0)
    assert _hex(fh.window_s) == _hex(W_B2)
    assert fh.sense is None
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = cli.main(["check", str(FHT_FILE)])
    printed = buf.getvalue()
    _p0(
        f"  (K) cli check rc (min, max) {_rc_all(rc, comm)}; family line present "
        f"{'family: f_human_birdcage' in printed}"
    )
    assert _rc_all(rc, comm) == (0, 0)
    # cli prints on rank 0 only
    assert comm.bcast("family: f_human_birdcage" in printed, root=0)

    # ---- (X) in-process one-key mutations -----------------------------------
    base = _mapping(FHT_FILE)

    def mut(fn):
        m = copy.deepcopy(base)
        fn(m)
        return m

    def _th(fn):
        def f(m):
            fn(m["thermal"])

        return f

    cases = [
        ("transient removed", _th(lambda t: t.pop("transient")), "thermal"),
        (
            "dt_s 0.0",
            _th(lambda t: t["transient"].__setitem__("dt_s", 0.0)),
            "thermal.transient.dt_s",
        ),
        (
            "n_steps 0",
            _th(lambda t: t["transient"].__setitem__("n_steps", 0)),
            "thermal.transient.n_steps",
        ),
        (
            "theta 1.5",
            _th(lambda t: t["transient"].__setitem__("theta", 1.5)),
            "thermal.transient.theta",
        ),
        ("balls removed", _th(lambda t: t.pop("balls")), "thermal"),
        (
            "mass_g 0.0",
            _th(lambda t: t["balls"].__setitem__("mass_g", 0.0)),
            "thermal.balls.mass_g",
        ),
        (
            "centres []",
            _th(lambda t: t["balls"].__setitem__("centres", [])),
            "thermal.balls.centres",
        ),
        ("sense up", _th(lambda t: t.__setitem__("sense", "up")), "thermal.sense"),
        (
            "unknown transient.dt",
            _th(lambda t: t["transient"].__setitem__("dt", 31.35)),
            "thermal.transient.dt",
        ),
    ]
    for label, fn, expected in cases:
        exc = _raises(mut(fn), f"<{label}>")
        _p0(f"  (X) {label}: field {exc.field!r} (expected {expected!r}); {exc.reason}")
        assert type(exc) is StudyConfigError, label
        assert exc.field == expected, (label, exc.field)
    m390 = _mapping(FH_FILE)
    m390["thermal"] = {"k_w_per_m_k": 0.6, "h_w_per_m2_k": 10.0, "t_inf_c": 20.0}
    exc = _raises(m390, "<390a mutation>")
    _p0(f"  (X) 390a's k/h/t_inf-only mutation: field {exc.field!r}; {exc.reason}")
    assert type(exc) is StudyConfigError and exc.field == "thermal"
    _p0(f"  (X) {len(cases) + 1} of {len(cases) + 1} mutations named their field")

    bad_dir = LOGS / "f_human_thermal_refused"
    _fresh_dir(bad_dir, comm)
    bad = bad_dir / "no_transient.yaml"
    if comm.rank == 0:
        bad_dir.mkdir(parents=True, exist_ok=True)
        bad.write_text(yaml.safe_dump(mut(_th(lambda t: t.pop("transient"))), sort_keys=False))
    comm.Barrier()
    out = bad_dir / "out"
    run_all = _rc_all(cli.main(["run", str(bad), "--output-dir", str(out)]), comm)
    check_all = _rc_all(cli.main(["check", str(bad)]), comm)
    exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
    _p0(
        f"  (X) cli on transient-less file: run rc {run_all}, check rc {check_all}, "
        f"out exists {exists}"
    )
    assert run_all == (2, 2) and check_all == (2, 2)
    assert not exists

    # ---- (c1) dropped key / (c2) one ULP ------------------------------------
    m1 = mut(_th(lambda t: t.pop("source_scale")))
    c1 = _load_from_mapping(m1).cooling
    default = RobinCooling(K_GEL, H_GEL, T_INF).source_scale
    d1 = [
        k
        for k in ("k_w_per_m_k", "h_w_per_m2_k", "t_inf_c", "source_scale", "degree")
        if _hex(getattr(c1, k)) != _hex(getattr(cool, k))
    ]
    _p0(
        f"  (c1) dropped source_scale diff {d1}; {_hex(c1.source_scale)} vs default "
        f"{_hex(default)}"
    )
    assert d1 == ["source_scale"], d1
    assert _hex(c1.source_scale) == _hex(default)
    G = float(js["G_per_w"])
    ulp = yaml.safe_load(repr(float(np.nextafter(G, 2.0 * G))))
    assert type(ulp) is float and ulp != G
    c2 = _load_from_mapping(mut(_th(lambda t: t.__setitem__("source_scale", ulp)))).cooling
    d2 = [
        k
        for k in ("k_w_per_m_k", "h_w_per_m2_k", "t_inf_c", "source_scale", "degree")
        if _hex(getattr(c2, k)) != _hex(getattr(cool, k))
    ]
    _p0(f"  (c2) one-ULP source_scale {ulp!r} diff {d2}")
    assert d2 == ["source_scale"], d2
    assert np.nextafter(G, 2.0 * G) == c2.source_scale


def _load_from_mapping(mapping):
    from fem_em_solver.study_yaml import _parse_mapping

    return _parse_mapping(mapping, "<mapping>")[3]


def test_wf1_fht_run_branch_no_mesh():
    """(Run-0): ``run_study_file``'s F-human branch with ``thermal:`` on fakes."""
    import yaml

    comm = MPI.COMM_WORLD
    s_node = np.array([[0.125 + 0.25j, -0.375 + 0.5j], [-0.375 + 0.5j, 0.625 - 0.75j]])
    y_in = {
        "ccw": np.array([0.01 + 0.02j, 0.03 - 0.04j]),
        "cw": np.array([-0.05 + 0.06j, 0.07 + 0.08j]),
    }
    gap, feed, vs = complex(0.2, -21.03), complex(7.6, -17.7), complex(1.02, 2.38)
    solves, temps = [], []

    def fake(study, comm_=MPI.COMM_WORLD):
        solves.append(1)
        study.result = SimpleNamespace(
            cells=369387,
            s_node=s_node,
            y_in=y_in,
            feeds=("P17", "P21"),
            az_feeds=(11.25, 101.25),
            az_mirror=56.25,
            gap_impedance_ohm=gap,
            feed_impedance_ohm=feed,
            source_voltage_v=vs,
        )
        return study.result

    def fake_temp(
        self, study, cooling, *, transient, balls, window_s, sense="ccw", comm=MPI.COMM_WORLD
    ):
        temps.append((cooling, transient, balls, window_s, sense))
        return SimpleNamespace(
            times_s=[0.0, 1.5, 3.0],
            volume_mean_k=[0.0, 0.25, 0.5],
            balls=[
                {
                    "centre": [0.01, 0.02, 0.0],
                    "radius_m": 0.0125,
                    "volume_m3": 8.0e-6,
                    "theta_k": [0.0, 0.375, 0.75],
                },
                {
                    "centre": [-0.01, 0.0, 0.05],
                    "radius_m": 0.0125,
                    "volume_m3": 7.5e-6,
                    "theta_k": [0.0, 0.125, 0.625],
                },
            ],
            six_minute={
                "window_s": 3.0,
                "volume_mean": {"peak_k": 0.25, "t_max_s": 3.0},
                "balls": [{"peak_k": 0.375, "t_max_s": 3.0}, {"peak_k": 0.21875, "t_max_s": 3.0}],
            },
            source_power_w=0.875,
            max_balance_residual=1.5e-10,
            cells=369387,
            n_sub=3306,
            sense=sense,
            cooling=cooling,
            transient=transient,
        )

    def boom(*args, **kwargs):
        raise AssertionError("no mesh / field solve in item 431a")

    run_dir = LOGS / "f_human_thermal_run0"
    _fresh_dir(run_dir, comm)
    with_th, without_th = run_dir / "with_thermal.yaml", run_dir / "without_thermal.yaml"
    if comm.rank == 0:
        run_dir.mkdir(parents=True, exist_ok=True)
        m = _mapping(FHT_FILE)
        with_th.write_text(yaml.safe_dump(m, sort_keys=False))
        del m["thermal"]
        without_th.write_text(yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(study_mod, "_study_solve_tuned_quadrature", fake)
        mp.setattr(TwoPortQuadratureDrive, "_temperature_port_study", fake_temp)
        mp.setattr(study_mod.Study, "_solve_port_fields", boom)
        mp.setattr(study_mod.FHumanBirdcageGeometry, "build", boom)
        run = run_study_file(with_th, output_dir=run_dir / "out", comm=comm)
        n_s = comm.allreduce(len(solves), op=MPI.MAX), comm.allreduce(len(solves), op=MPI.MIN)
        n_t = comm.allreduce(len(temps), op=MPI.MAX), comm.allreduce(len(temps), op=MPI.MIN)
        solves.clear()
        temps_before = len(temps)
        run2 = run_study_file(without_th, output_dir=run_dir / "out2", comm=comm)
        n_s2 = comm.allreduce(len(solves), op=MPI.MAX), comm.allreduce(len(solves), op=MPI.MIN)
        n_t2 = comm.allreduce(len(temps) - temps_before, op=MPI.MAX)
    summary = run.summary
    _p0(
        f"\n[WF-1 s10a] (Run-0) solve calls (max, min) {n_s}; temperature calls {n_t}; "
        f"keys {sorted(summary)}"
    )
    assert n_s == (1, 1) and n_t == (1, 1)
    cooling, transient, balls, window_s, sense = temps[0]
    js = _series_json()
    assert js is not None
    assert _hex(cooling.source_scale) == _hex(js["G_per_w"])
    assert _hex(transient.dt_s) == _hex(T_STAR / N_MAIN) and transient.n_steps == N_MAIN
    assert [_hex(v) for v in balls["centres"][0]] == [_hex(v) for v in js["c_star_m"]]
    assert _hex(window_s) == _hex(W_B2) and sense == "ccw"
    exp_t = {
        "sense": "ccw",
        "times_s": [0.0, 1.5, 3.0],
        "volume_mean_k": [0.0, 0.25, 0.5],
        "balls": [
            {
                "centre": [0.01, 0.02, 0.0],
                "radius_m": 0.0125,
                "volume_m3": 8.0e-6,
                "theta_k": [0.0, 0.375, 0.75],
            },
            {
                "centre": [-0.01, 0.0, 0.05],
                "radius_m": 0.0125,
                "volume_m3": 7.5e-6,
                "theta_k": [0.0, 0.125, 0.625],
            },
        ],
        "six_minute": {
            "window_s": 3.0,
            "volume_mean": {"peak_k": 0.25, "t_max_s": 3.0},
            "balls": [{"peak_k": 0.375, "t_max_s": 3.0}, {"peak_k": 0.21875, "t_max_s": 3.0}],
        },
        "source_power_w": 0.875,
        "max_balance_residual": 1.5e-10,
        "source_scale": float(js["G_per_w"]),
    }
    _p0(f"  (Run-0) summary temperature == fake: {summary.get('temperature') == exp_t}")
    assert summary["temperature"] == exp_t
    assert summary["family"] == "f_human_birdcage" and summary["cells"] == 369387
    assert run.thermal is not None and run.thermal.sense == "ccw"
    same = True
    if comm.rank == 0:
        with open(run.paths["summary"]) as fh:
            same = json.load(fh) == summary
    same = comm.allreduce(bool(same), op=MPI.LAND)
    _p0(f"  (Run-0) written JSON == in-memory summary: {same}")
    assert same
    _p0(
        f"  (Run-0) without thermal: solve calls {n_s2}, temperature calls (max) {n_t2}, "
        f"'temperature' in summary {'temperature' in run2.summary}"
    )
    assert n_s2 == (1, 1) and n_t2 == 0
    assert "temperature" not in run2.summary and run2.thermal is None
