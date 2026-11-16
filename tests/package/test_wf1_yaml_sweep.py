"""`WF-1` step 21: ``study.sweep_hz`` on item 192's four-port F-small file drives
the `WF-16` ``Study.solve_sweep`` through ``femem run`` — one mesh, three
frequencies, one ``.s4p``, each S equal to the `ANS-4` exact-strip ``P1`` record.

Brief: docs/planning/chunks/WF-1.md, "§9 item 572" and the review ruling
2026-10-20 10:30.  The solve legs run in the complex build and only with
``FEM_EM_WF1_S21=1`` (then they assert the complex build and the record's
width, never skip).  The parse node builds no mesh and runs in either build.

Anchors: (K) the file's Study = ``f_small_four_port.yaml``'s, the parsed sweep
= ``FREQUENCIES`` by ``float.hex``; (X) refusals naming their field (shape,
entries, order, conductor, hole, capacitor sheets, sphere, F-human, SAR / B₁⁺ /
probes / tune, the Palace backend); (R′) per-frequency S vs ``P1`` at
``LEG_D0_REPRODUCTION_BAND`` on ``ANS4_EXACT_STRIP_CELL_RECORD`` cells, one mesh;
(T′) the one Touchstone vs memory at ``ANCHOR_TOL`` / ``FREQ_ATOL_HZ``; (G′)
reciprocity and passivity per frequency; (U′) JSON = memory, no scalar solve;
(I-Py) θ_file = θ_Python at ``L0_BAND``.  Control (rule (e), asserted): (c-swap)
the 10 ↔ 128 MHz exchange misses.  Every band, record and helper is imported.
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
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study_yaml import _load, _run_study_file_palace, _study_from_mapping

from tests.complex_mode import IS_COMPLEX
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ
from tests.package.test_wf16_study_four_port import _rel
from tests.package.test_wf16_study_sweep import C4_HAND_TABLE, KEYS, _port_study, _records
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    FREQUENCIES,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

S21_ENV = "FEM_EM_WF1_S21"
SWEEP = WF1 / "f_small_sweep.yaml"
FOUR_PORT = WF1 / "f_small_four_port.yaml"
NPZ = LOGS / "s21" / "sweep-n2.npz"


def _knob() -> bool:
    raw = os.environ.get(S21_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S21_ENV}={raw!r}: only unset/empty or '1' are accepted")


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


def _with_sweep(mapping, sweep):
    return {**mapping, "study": {**mapping["study"], "sweep_hz": sweep}}


GOOD_SWEEP = [10000000.0, 64000000.0, 128000000.0]
BAD_LISTS = (
    ([], "study.sweep_hz"),
    (6.4e7, "study.sweep_hz"),
    ([64000000.0, 10000000.0], "study.sweep_hz"),
    ([10000000.0, 10000000.0], "study.sweep_hz"),
    ([0.0], "study.sweep_hz[0]"),
    ([-1.0], "study.sweep_hz[0]"),
    (["10 MHz"], "study.sweep_hz[0]"),
)


# ------------------------------------------------------------------ no mesh


def test_wf1_s21_parse_and_refusals(monkeypatch, capsys):
    """(K), (X): no mesh, no solve, either build."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                ThermalBallGeometry, ThermalBoxGeometry):
        monkeypatch.setattr(cls, "build", stub)

    diff = _spec_diff(load_study(SWEEP), load_study(FOUR_PORT))
    study, kind, outputs, cooling = _load(SWEEP)
    got = [float(f).hex() for f in outputs["sweep_hz"]]
    want = [float(f).hex() for f in FREQUENCIES.values()]
    _say(capsys, f"\n[WF-1 s21] (K) spec diff {diff}; kind {kind!r}; outputs keys "
         f"{sorted(outputs)}; sweep_hz {outputs['sweep_hz']}; hex equal {got == want}; "
         f"frequency_hz {study.frequency_hz!r}")
    assert diff == [], diff
    assert kind == "four_port_birdcage"
    assert cooling is None
    assert set(outputs) == {"sweep_hz"}
    assert isinstance(outputs["sweep_hz"], tuple)
    assert got == want, (got, want)

    capsys.readouterr()
    rc = cli.main(["check", str(SWEEP)])
    lines = capsys.readouterr().out.splitlines()
    rc_all = _rc_all(rc, comm)
    _say(capsys, f"  (K) femem check {SWEEP.name}: rc {rc_all}")
    assert rc_all == (0, 0)
    if comm.rank == 0:
        assert [ln for ln in lines if ln.startswith("family: ")] == [
            "family: four_port_birdcage"
        ]

    base = _mapping(SWEEP)
    sar_b1 = _mapping(WF1 / "f_small_sar_b1.yaml")["outputs"]
    probes = _mapping(WF1 / "f_small_em_probes.yaml")["outputs"]["probes"]
    rows = [(_with_sweep(base, s), f) for s, f in BAD_LISTS]
    for name in ("f_small_copper.yaml", "f_small_hole_pec_control.yaml",
                 "f_small_hole_pec_capacitors.yaml", "sphere_mat4.yaml",
                 "f_human_hole_copper_tune.yaml"):
        rows.append((_with_sweep(_mapping(WF1 / name), GOOD_SWEEP), "study.sweep_hz"))
    rows += [
        ({**base, "outputs": {"sar": sar_b1["sar"]}}, "study.sweep_hz"),
        ({**base, "outputs": {"b1_plus": sar_b1["b1_plus"]}}, "study.sweep_hz"),
        ({**base, "outputs": {"probes": probes}}, "study.sweep_hz"),
        ({**base, "tune": {}}, "study.sweep_hz"),
        (base, None),
    ]
    fields = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for i, ((_, want_f), g) in enumerate(zip(rows, fields)):
        _say(capsys, f"  (X) row {i}: want {want_f!r} got {g!r}")
    assert fields == [f for _, f in rows], [
        (f, g) for (_, f), g in zip(rows, fields) if f != g
    ]

    out = LOGS / "s21" / "palace_refused"
    _fresh_dir(out, comm)
    g_pal = _field(lambda: _run_study_file_palace(SWEEP, output_dir=out, comm=comm))
    written = _all(comm, not out.exists())
    _say(capsys, f"  (X) --backend palace on {SWEEP.name}: field {g_pal!r}; nothing written "
         f"{written}")
    assert g_pal == "--backend"
    assert written

    py, _ = _port_study(1)
    parity = []
    for s, _f in BAD_LISTS:
        if not isinstance(s, list) or any(isinstance(v, str) for v in s):
            continue  # not a list of numbers: a YAML-only shape
        with pytest.raises(ValueError):
            py.solve_sweep(tuple(s), comm)
        parity.append(s)
    _say(capsys, f"  (X) Python solve_sweep raises ValueError on {parity}")
    assert len(parity) == 5  # [], decreasing, duplicate, [0.0], [-1.0]
    _say(capsys, f"  build calls {calls['n']}")
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _require():
    if not _knob():
        pytest.skip(f"{S21_ENV} unset: the sweep solve legs are environment-selected "
                    "(complex build)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port port Study needs the complex build"
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        comm.size, ANS4_EXACT_STRIP_RECORD_RANK_WIDTH)
    return comm


def _unpairs(m):
    return np.array([[complex(re, im) for re, im in row] for row in m], dtype=np.complex128)


def test_wf1_s21_sweep_file_n2():
    comm = _require()
    t00 = time.perf_counter()
    out = LOGS / "s21" / "sweep_n2"
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    run = run_study_file(SWEEP, output_dir=out, comm=comm)
    t_run = time.perf_counter() - t0
    s = run.summary
    sw = run.result
    _p0(f"\n[WF-1 s21] -n {comm.size}: run_study_file {t_run:.1f} s (summary elapsed "
        f"{s['elapsed_s']:.1f} s); mesh {s['mesh_time_s']:.2f} s; sweeps "
        + " / ".join(f"{k} MHz {t:.2f} s" for k, t in zip(KEYS, s["sweep_times_s"]))
        + f"; {s['cells']} cells; frequency_hz {s['frequency_hz']!r} beside sweep_hz "
        f"{s['sweep_hz']}")

    # (R′) per frequency, one mesh
    rec = _records(1)
    s_file = np.array([_unpairs(m) for m in s["s_matrices"]])
    rels = {k: _rel(s_file[i], rec[k]) for i, k in enumerate(KEYS)}
    one_mesh = all(r.mesh is sw.results[0].mesh for r in sw.results)
    _p0("  (R') ||S_f - P1_f||/||P1_f|| " + " / ".join(f"{k} MHz {v:.3e}" for k, v in rels.items())
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}); cells {s['cells']} (record "
        f"{ANS4_EXACT_STRIP_CELL_RECORD}); one mesh {one_mesh}")
    assert s["cells"] == ANS4_EXACT_STRIP_CELL_RECORD
    assert one_mesh
    assert len(sw.results) == len(KEYS)
    assert [float(f).hex() for f in s["sweep_hz"]] == [float(f).hex() for f in FREQUENCIES.values()]
    for k, v in rels.items():
        assert v <= LEG_D0_REPRODUCTION_BAND, k

    # (c-swap) asserted: the 10 <-> 128 MHz exchange misses
    i10, i128 = KEYS.index("10"), KEYS.index("128")
    sw_10 = _rel(s_file[i10], rec["128"])
    sw_128 = _rel(s_file[i128], rec["10"])
    _p0(f"  (c-swap) S_10 vs P1_128 {sw_10:.4f} (hand {C4_HAND_TABLE[1][('10', '128')]}); "
        f"S_128 vs P1_10 {sw_128:.4f} (hand {C4_HAND_TABLE[1][('128', '10')]}) "
        f"(ASSERTED > {LEG_D0_REPRODUCTION_BAND:.0e})")
    assert sw_10 > LEG_D0_REPRODUCTION_BAND
    assert sw_128 > LEG_D0_REPRODUCTION_BAND

    # (T′) the one Touchstone
    ts = run.paths["touchstone"]
    freqs, mats, ports = load_touchstone(ts)
    loaded = np.asarray(mats, dtype=np.complex128)
    d_ts = [float(np.max(np.abs(loaded[i] - sw.s_matrices[i]))) for i in range(len(KEYS))]
    d_f = [abs(float(freqs[i]) - f) for i, f in enumerate(s["sweep_hz"])]
    _p0(f"  (T') {ts.name}: {len(freqs)} points; |df| {d_f}; max|dS| "
        + " / ".join(f"{v:.3e}" for v in d_ts) + f" (ASSERTED <= {ANCHOR_TOL:.0e})")
    assert ts.name == f"{SWEEP.stem}.s4p" and s["touchstone"] == ts.name
    assert len(freqs) == len(KEYS) and loaded.shape == (len(KEYS), 4, 4)
    assert all(d <= FREQ_ATOL_HZ for d in d_f)
    assert all(d <= ANCHOR_TOL for d in d_ts)

    # (G′) PORT-11's identities per frequency
    for i, k in enumerate(KEYS):
        S = sw.s_matrices[i]
        recip = float(np.linalg.norm(S - S.T) / np.linalg.norm(S))
        smax = float(np.linalg.svd(S, compute_uv=False).max())
        _p0(f"  (G') {k} MHz reciprocity {recip:.3e} (ASSERTED <= {RECIPROCITY_BAND:.0e}); "
            f"sigma_max {smax:.6f} (ASSERTED <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e})")
        assert recip <= RECIPROCITY_BAND, k
        assert smax <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, k

    # (U′) JSON = memory; summary = the sweep exactly; no scalar solve
    # Timings are rank-local (each rank's own clock); rank 0's are the ones on disk.
    # First W1a (20261020T201341Z) compared them on rank 1 too and read False.
    timing = ("elapsed_s", "mesh_time_s", "sweep_times_s")
    js = json.loads(run.paths["summary"].read_text())
    rest = {k: v for k, v in s.items() if k not in timing}
    js_rest = {k: v for k, v in js.items() if k not in timing}
    diff_keys = sorted(k for k in set(js) | set(s) if js.get(k) != s.get(k))
    _p0(f"  (U') rank 0 keys differing JSON vs memory: {diff_keys}")
    same_json = _all(comm, js_rest == rest and (comm.rank != 0 or js == s))
    exact = bool(np.array_equal(s_file, sw.s_matrices))
    _p0(f"  (U') JSON == memory {same_json}; summary S == sweep S exactly {exact}; "
        f"study.result is None {run.study.result is None}")
    assert same_json
    assert exact
    assert run.study.result is None

    agree = all(np.array_equal(a, s_file) for a in comm.allgather(s_file))
    assert agree
    if comm.rank == 0:
        NPZ.parent.mkdir(parents=True, exist_ok=True)
        np.savez(NPZ, frequencies=np.asarray(s["sweep_hz"]), s=s_file,
                 cells=np.array(s["cells"]))
    comm.Barrier()
    _p0(f"  wrote {NPZ}; node {time.perf_counter() - t00:.1f} s")


def test_wf1_s21_sweep_python_twin():
    comm = _require()
    exists = comm.allreduce(bool(NPZ.exists()), op=MPI.LAND)
    assert exists, f"{NPZ} absent: run test_wf1_s21_sweep_file_n2 at -n 2 first (never re-created here)"
    ref = np.load(NPZ)
    t00 = time.perf_counter()
    study, _ = _port_study(1)
    sw = study.solve_sweep(tuple(FREQUENCIES.values()), comm)
    _p0(f"\n[WF-1 s21] (I-Py) -n {comm.size}: Python solve_sweep {time.perf_counter() - t00:.1f} s; "
        f"mesh {sw.mesh_time:.2f} s; sweeps "
        + " / ".join(f"{k} MHz {t:.2f} s" for k, t in zip(KEYS, sw.sweep_times))
        + f"; {sw.results[0].cells} cells (file {int(ref['cells'])})")
    assert np.array_equal(np.asarray(sw.frequencies_hz), ref["frequencies"])
    assert int(sw.results[0].cells) == int(ref["cells"])
    rels = [_rel(sw.s_matrices[i], ref["s"][i]) for i in range(len(KEYS))]
    _p0("  (I-Py) ||S_py - S_file||/||S_file|| " + " / ".join(
        f"{k} MHz {v:.3e}" for k, v in zip(KEYS, rels)) + f" (ASSERTED <= {L0_BAND:.0e})")
    for k, v in zip(KEYS, rels):
        assert v <= L0_BAND, k
