"""`WF-15` step 11 (§9 item 728, §10 U109, provisional under rule (6)): the rod ``Study``
on a sweep — ``FourPortBirdcageGeometry(implant_rod=ROD)``, degree 1, ``LumpedPortDrive``.

Knob ``FEM_EM_WF15_S11`` = ``1``; unset, the coding test runs (no mesh) and the solve node
skips.  Every band is imported; 128 MHz has no record, so its σ_max is printed.

Asserted: ``Study.solve_sweep`` at (64, 128) MHz equals one-frequency rod ``Study.solve``
at each f, per entry, at ``IDENTITY_BAND``; (i) ``RECIPROCITY_BAND`` and (ii) σ_max ≤ 1 +
``PASSIVITY_SIGMA_TOLERANCE``, column power ≤ 1, at each f; cells 168 507 on every build;
the 64 MHz σ_max reproduces item 709's print to half its last digit.  Control: the sweep's
64 MHz matrix misses the 128 MHz one-frequency solve by > ``IDENTITY_BAND``, so a sweep
that returned one frequency twice fails.  YAML: ``implant:`` beside ``study.sweep_hz``
parses, with no second solve.

Not claimed: a tip Q at 128 MHz, a resonance, the two-port drive, a tissue model.
"""

from __future__ import annotations

import dataclasses
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.study_yaml import _load as _yaml_load

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf15_study_implant_rod import (
    RECORD_WIDTH,
    ROD_EXACT_CELL_RECORD,
    _mapping,
    _phantom,
    _rod_yaml,
    _study,
)
from tests.validation.test_imp3_rod_birdcage import _matrix_quantities
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ, FREQUENCY_128_HZ
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KNOB = "FEM_EM_WF15_S11"
FREQS = (FREQUENCY_64_HZ, FREQUENCY_128_HZ)
#: Item 709's σ_max at 64 MHz, `20261101T104952Z_WF-15-s10-W1.log`, printed at `.12f`.
ITEM_709_SIGMA = 0.999767876817
SIGMA_PRINT_BOUND = 5e-13  # half the last printed digit

needs_knob = pytest.mark.skipif(os.environ.get(KNOB) != "1", reason=f"{KNOB} unset")
coding = pytest.mark.skipif(os.environ.get(KNOB) == "1", reason=f"{KNOB} set: coding tests skip")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


@coding
def test_s11_yaml_rod_sweep_parses(tmp_path):
    """`implant:` and `study.sweep_hz` together parse, with no mesh and no solve."""
    import yaml

    comm = MPI.COMM_WORLD
    base = _mapping(implant=_rod_yaml())
    mapping = {**base, "study": {**base["study"], "sweep_hz": [float(f) for f in FREQS]}}
    path = tmp_path / f"rod-sweep-{comm.rank}.yaml"
    path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    study, kind, outputs, _cooling = _yaml_load(path)
    assert kind == "four_port_birdcage"
    assert study.geometry == FourPortBirdcageGeometry(implant_rod=ROD)
    assert study.geometry.implant_rod is not None
    assert [float(f).hex() for f in outputs["sweep_hz"]] == [float(f).hex() for f in FREQS]
    assert set(outputs) == {"sweep_hz"}, sorted(outputs)
    _p0(f"\n[WF-15 s11] YAML admits the rod sweep: sweep_hz {outputs['sweep_hz']}, "
        f"implant_rod keys {sorted(study.geometry.implant_rod)}")


@pytest.fixture(scope="module")
def s11():
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_WIDTH, f"the record holds at -n {RECORD_WIDTH}, got -n {comm.size}"
    geom = FourPortBirdcageGeometry(implant_rod=ROD)
    t0 = time.perf_counter()
    sweep = _study(geom).solve_sweep(FREQS, comm)
    t_sweep = time.perf_counter() - t0
    t0 = time.perf_counter()
    points = {}
    for f in FREQS:
        material = dataclasses.replace(_phantom(), frequency_hz=f)
        points[f] = Study(geom, material, LumpedPortDrive(), frequency_hz=f, degree=1).solve(comm)
    t_points = time.perf_counter() - t0
    _p0(f"\n[WF-15 s11] -n {comm.size}: sweep {t_sweep:.1f} s (mesh {sweep.mesh_time:.1f}, solves "
        f"{[round(t, 1) for t in sweep.sweep_times]}); two point Studies {t_points:.1f} s")
    return {"sweep": sweep, "points": points}


@complex_only
@needs_knob
def test_s11_all_gates(s11):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    sweep, points = s11["sweep"], s11["points"]
    gates = {}
    for f, res in zip(FREQS, sweep.results):
        mhz = f / 1e6
        point = points[f]
        mq = _matrix_quantities(res)
        sigma = float(mq["sigma"][0])
        rel = _rel(res.s_matrix, point.s_matrix)
        entry = float(np.max(np.abs(np.asarray(res.s_matrix) - np.asarray(point.s_matrix))
                             / np.abs(np.asarray(point.s_matrix))))
        gates[f"{mhz:g} S"] = rel <= IDENTITY_BAND and entry <= IDENTITY_BAND
        gates[f"{mhz:g} i"] = mq["recip"] <= RECIPROCITY_BAND
        gates[f"{mhz:g} ii_sigma"] = sigma <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
        gates[f"{mhz:g} ii_colpow"] = bool(np.all(mq["colpow"] <= 1.0))
        gates[f"{mhz:g} cells sweep"] = res.cells == ROD_EXACT_CELL_RECORD
        gates[f"{mhz:g} cells point"] = point.cells == ROD_EXACT_CELL_RECORD
        _p0(f"[WF-15 s11] {mhz:g} MHz: S sweep vs point rel {rel:.3e}, worst entry {entry:.3e} "
            f"(IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED); (i) {mq['recip']:.3e} (band "
            f"{RECIPROCITY_BAND:.0e}); (ii) sigma_max {sigma:.12f}, colpow "
            f"{np.array2string(mq['colpow'], precision=7)}; cells {res.cells} / {point.cells}")
    sig64 = float(_matrix_quantities(sweep.results[0])["sigma"][0])
    gates["sigma 64 record"] = abs(sig64 - ITEM_709_SIGMA) <= SIGMA_PRINT_BOUND
    _p0(f"[WF-15 s11] 64 MHz sigma_max vs item 709 {ITEM_709_SIGMA:.12f}: |diff| "
        f"{abs(sig64 - ITEM_709_SIGMA):.3e} (bound {SIGMA_PRINT_BOUND:g}, ASSERTED)")
    miss = _rel(sweep.results[0].s_matrix, points[FREQS[1]].s_matrix)
    gates["control"] = miss > IDENTITY_BAND
    _p0(f"[WF-15 s11] control: sweep 64 MHz vs 128 MHz point miss {miss:.4e} (asserted > "
        f"{IDENTITY_BAND:.0e}), {miss / IDENTITY_BAND:.3e} x the band (POOLED_SEPARATION_FLOOR "
        f"{POOLED_SEPARATION_FLOOR:.0f}x, printed)")
    red = [k for k, v in gates.items() if not v]
    _p0(f"[WF-15 s11] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red
