"""`WF-16` step 18 (§9 item 780, §10 U115): ``Study.solve_sweep(..., fields_at=)``.

One keyword on the existing sweep: ``fields_at`` names one frequency of
``frequencies_hz`` whose ``results[i].fields`` keeps the drive fields, taken
inside that frequency's own ``run_n_port_sparameter_sweep`` call
(``keep_fields=True``); every other frequency keeps ``fields=None``.

Fixture: `mri:38`'s rod four-port (``FourPortBirdcageGeometry(implant_rod=ROD)``,
degree 1, ``LumpedPortDrive``), swept at (64, 128) MHz with ``fields_at=64 MHz``.
Knob ``FEM_EM_WF16_S18`` = ``1``; unset, only the no-mesh refusal test runs.

Asserted (collected first, then asserted):

* (I) ``|B1+|`` (``Study.b1_plus_at``, every drive, one point batch) and the 10 g
  SAR (``Study.sar``, drive P1, one centre) from the ``fields_at`` result equal
  those of the private ``Study._solve_port_fields`` route (a separate solve on
  its own mesh build) within ``IDENTITY_BAND``.  Both sides go through the same
  public post-processing: the swept result is set as a fresh Study's
  ``fields_result``, so the operator code is shared and only the solve differs.
  Points are sampled in the same batch on both sides (``.eval`` bytes depend on
  batch composition).  A red is read first against `OPS-77` / `OPS-78`'s
  measured separate-solve drift (largest 1.147e-13).
* (S0) one ``LinearProblem`` construction and one solve per frequency (the
  `PORT-19` reused factorisation), counted by a ``_CountingLinearProblem`` spy
  (replicated from ``tests/validation/test_port19_factor_reuse.py:132-181``),
  equal per frequency to the ``fields_at=None`` run's counts.
* (R) a ``fields_at`` not in ``frequencies_hz`` (float equality; also NaN and a
  non-number) raises ``ValueError`` naming ``fields_at`` before any mesh.
* Control: ``fields_at=None`` leaves every ``fields`` as ``None``; with
  ``fields_at`` set, only the 64 MHz result holds fields and only its sweep
  call received ``keep_fields=True``.

Not claimed: fields at every frequency, conductor sweeps (still refused).
"""

from __future__ import annotations

import dataclasses
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.core.time_harmonic as th_module
import fem_em_solver.ports.sparameters as sparam_module
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf15_study_implant_rod import RECORD_WIDTH, _phantom, _study
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ, FREQUENCY_128_HZ

KNOB = "FEM_EM_WF16_S18"
FREQS = (FREQUENCY_64_HZ, FREQUENCY_128_HZ)
FIELDS_AT = FREQUENCY_64_HZ
#: `PORT-19` step 2: the reused lumped-sheet sweep constructs and solves one
#: ``LinearProblem`` per frequency (``EXPECTED_REUSED_FACTORISATIONS`` there).
EXPECTED_FACTORISATIONS_PER_FREQUENCY = 1
#: `OPS-77`'s largest separate-solve sample (`20261102T122110Z_OPS-77-reader.log:36-44`),
#: printed beside an (I) reading, never asserted.
OPS77_LARGEST_DRIFT = 1.147e-13
#: Points inside the phantom, away from the rod (centre (0.015, 0, 0), radius 1 mm).
POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [-0.01, 0.0, 0.0],
        [0.0, 0.01, 0.0],
        [0.0, -0.01, 0.01],
        [0.005, 0.005, -0.01],
    ],
    dtype=np.float64,
)
SAR_CENTRE = (0.0, 0.0, 0.0)

needs_knob = pytest.mark.skipif(os.environ.get(KNOB) != "1", reason=f"{KNOB} unset")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def test_wf16_s18_fields_at_refused_before_mesh(monkeypatch):
    """(R) ``fields_at`` outside ``frequencies_hz`` raises ``ValueError`` before any mesh."""

    def _no_build(self, comm):
        raise AssertionError("build reached: the refusal must come first")

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _no_build)
    study = _study(FourPortBirdcageGeometry(implant_rod=ROD))
    bad = (100e6, FREQUENCY_64_HZ * (1 + 1e-15), float("nan"), "64 MHz", 0.0)
    for value in bad:
        with pytest.raises(ValueError, match="fields_at"):
            study.solve_sweep(FREQS, MPI.COMM_WORLD, fields_at=value)
    assert getattr(study, "sweep_result", None) is None
    # The in-set value passes the check and reaches the build (the patched raise).
    with pytest.raises(AssertionError, match="build reached"):
        study.solve_sweep(FREQS, MPI.COMM_WORLD, fields_at=FIELDS_AT)
    _p0(
        f"\n[WF-16 s18] (R) {len(bad)} values outside {FREQS} refused naming fields_at, "
        "before any mesh; the in-set value reaches build"
    )


class _Counters:
    def __init__(self):
        self.init = 0
        self.solve = 0
        self.calls = []  # (init, solve, keep_fields) per sweep call


def _spied_sweep(study, comm, fields_at):
    """``study.solve_sweep`` with the factorisation spy, counts per frequency."""
    counters = _Counters()
    original_lp = th_module.LinearProblem
    original_sweep = sparam_module.run_n_port_sparameter_sweep

    class _CountingLinearProblem(original_lp):
        def __init__(self, *args, **kwargs):
            counters.init += 1
            super().__init__(*args, **kwargs)

        def solve(self):
            counters.solve += 1
            return super().solve()

    def _counted_sweep(*args, **kwargs):
        counters.init = counters.solve = 0
        out = original_sweep(*args, **kwargs)
        counters.calls.append(
            (counters.init, counters.solve, bool(kwargs.get("keep_fields", False)))
        )
        return out

    th_module.LinearProblem = _CountingLinearProblem
    sparam_module.run_n_port_sparameter_sweep = _counted_sweep
    try:
        comm.Barrier()
        t0 = time.perf_counter()
        sweep = study.solve_sweep(FREQS, comm, fields_at=fields_at)
        comm.Barrier()
        elapsed = time.perf_counter() - t0
    finally:
        th_module.LinearProblem = original_lp
        sparam_module.run_n_port_sparameter_sweep = original_sweep
    reduced = []
    for init, solve, keep in counters.calls:
        reduced.append(
            (
                (int(comm.allreduce(init, op=MPI.MIN)), int(comm.allreduce(init, op=MPI.MAX))),
                (int(comm.allreduce(solve, op=MPI.MIN)), int(comm.allreduce(solve, op=MPI.MAX))),
                keep,
            )
        )
    return sweep, reduced, elapsed


def _readings(study, comm):
    """``|B1+|`` (every drive, one batch) and 10 g SAR (P1, one centre) via the public methods."""
    b1 = study.b1_plus_at(POINTS, comm=comm)
    sar = study.sar(mass_g=10.0, centres=[SAR_CENTRE], drives=["P1"], comm=comm)
    return b1, sar


def _rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


@pytest.fixture(scope="module")
def s18():
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_WIDTH, f"the record holds at -n {RECORD_WIDTH}, got -n {comm.size}"
    geom = FourPortBirdcageGeometry(implant_rod=ROD)
    keep, keep_counts, t_keep = _spied_sweep(_study(geom), comm, FIELDS_AT)
    none, none_counts, t_none = _spied_sweep(_study(geom), comm, None)
    _p0(
        f"\n[WF-16 s18] -n {comm.size}: sweep fields_at={FIELDS_AT:g} {t_keep:.1f} s (solves "
        f"{[round(t, 1) for t in keep.sweep_times]}); sweep fields_at=None {t_none:.1f} s"
    )
    t0 = time.perf_counter()
    route = Study(
        geom,
        dataclasses.replace(_phantom(), frequency_hz=FIELDS_AT),
        LumpedPortDrive(),
        frequency_hz=FIELDS_AT,
        degree=1,
    )
    route_res = route._solve_port_fields(comm)
    _p0(
        f"[WF-16 s18] private _solve_port_fields route {time.perf_counter() - t0:.1f} s "
        f"(cells {route_res.cells})"
    )
    # The swept result through the same public post-processing: a fresh Study
    # whose cached field-keeping result is the sweep's fields_at result.
    swept = _study(geom)
    swept.fields_result = keep.results[FREQS.index(FIELDS_AT)]
    t0 = time.perf_counter()
    read_route = _readings(route, comm)
    read_swept = _readings(swept, comm)
    _p0(
        f"[WF-16 s18] readings (both sides) {time.perf_counter() - t0:.1f} s; omega route "
        f"{route.omega!r} vs swept {swept.omega!r}"
    )
    return {
        "keep": keep,
        "none": none,
        "keep_counts": keep_counts,
        "none_counts": none_counts,
        "route_res": route_res,
        "read_route": read_route,
        "read_swept": read_swept,
    }


@complex_only
@needs_knob
def test_wf16_s18_fields_at_gates(s18):
    """(I), (S0) and the control, collected before asserting so one red hides no other."""
    gates = {}
    keep, none = s18["keep"], s18["none"]
    # Control: the old contract, and only the named frequency holds fields.
    gates["control none"] = all(r.fields is None for r in none.results)
    held = [r.fields is not None for r in keep.results]
    gates["control fields_at"] = held == [f == FIELDS_AT for f in FREQS]
    _p0(
        f"[WF-16 s18] control: fields_at=None fields {[r.fields for r in none.results]}; "
        f"fields_at={FIELDS_AT:g} held {held} for {FREQS}; fields ids "
        f"{sorted(keep.results[0].fields or {})}"
    )
    # (S0)
    kc, nc = s18["keep_counts"], s18["none_counts"]
    one = (EXPECTED_FACTORISATIONS_PER_FREQUENCY, EXPECTED_FACTORISATIONS_PER_FREQUENCY)
    gates["S0 one per frequency"] = len(kc) == len(FREQS) and all(
        c[0] == one and c[1] == one for c in kc
    )
    gates["S0 equals None run"] = [c[:2] for c in kc] == [c[:2] for c in nc]
    gates["S0 keep only at fields_at"] = [c[2] for c in kc] == [
        f == FIELDS_AT for f in FREQS
    ] and not any(c[2] for c in nc)
    _p0(
        f"[WF-16 s18] (S0) per frequency (init (min,max), solve (min,max), keep_fields): "
        f"fields_at {kc}; None {nc}; expected {EXPECTED_FACTORISATIONS_PER_FREQUENCY} (ASSERTED)"
    )
    # (I)
    b1_r, sar_r = s18["read_route"]
    b1_s, sar_s = s18["read_swept"]
    gates["I drives"] = sorted(b1_r) == sorted(b1_s)
    for pid in sorted(b1_r):
        vr, okr = b1_r[pid]
        vs, oks = b1_s[pid]
        n_valid = int(np.count_nonzero(okr))
        gates[f"I B1 {pid} valid"] = bool(np.array_equal(okr, oks)) and n_valid >= 1
        rel = _rel(vs[okr], vr[okr]) if n_valid else float("inf")
        gates[f"I B1 {pid}"] = rel <= IDENTITY_BAND
        _p0(
            f"[WF-16 s18] (I) |B1+| {pid}: valid {n_valid}/{len(POINTS)}, max rel {rel:.3e} "
            f"(IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED; OPS-77 largest drift "
            f"{OPS77_LARGEST_DRIFT:.3e}, printed); route {np.array2string(vr, precision=9)}"
        )
    sr = float(sar_r["P1"][0]["averaged_sar_w_per_kg"])
    ss = float(sar_s["P1"][0]["averaged_sar_w_per_kg"])
    rel_sar = abs(ss - sr) / abs(sr)
    gates["I SAR P1"] = rel_sar <= IDENTITY_BAND and sr > 0.0
    _p0(
        f"[WF-16 s18] (I) SAR10g P1 at {SAR_CENTRE}: route {sr:.12e} W/kg, swept {ss:.12e} W/kg, "
        f"rel {rel_sar:.3e} (IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED)"
    )
    s_rel = float(
        np.linalg.norm(np.asarray(keep.results[0].s_matrix) - s18["route_res"].s_matrix)
        / np.linalg.norm(s18["route_res"].s_matrix)
    )
    s_none = float(
        np.linalg.norm(np.asarray(keep.results[0].s_matrix) - none.results[0].s_matrix)
        / np.linalg.norm(none.results[0].s_matrix)
    )
    _p0(
        f"[WF-16 s18] printed: 64 MHz S swept vs route rel {s_rel:.3e}; fields_at vs None sweep "
        f"rel {s_none:.3e}"
    )
    red = sorted(k for k, v in gates.items() if not v)
    red_all = sorted(set(MPI.COMM_WORLD.allreduce(red, op=MPI.SUM)))
    _p0(f"[WF-16 s18] gates red: {red_all or 'none'} of {len(gates)}")
    assert not red_all, red_all
