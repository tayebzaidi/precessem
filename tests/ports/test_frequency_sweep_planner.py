"""Tests for frequency sweep planning utility (chunk D4)."""

from __future__ import annotations

import itertools
import math
import random
from decimal import Decimal
from fractions import Fraction
from typing import Sequence

import numpy as np
import pytest

from fem_em_solver.ports import FrequencySweepPlan, plan_frequency_sweep


def test_plan_frequency_sweep_coarse_only_is_deterministic_and_inclusive():
    plan = plan_frequency_sweep(
        start_hz=120.0e6,
        stop_hz=130.0e6,
        coarse_step_hz=5.0e6,
        z0_ohm=50.0,
        port_order=("P1", "P2", "P3", "P4"),
    )

    assert isinstance(plan, FrequencySweepPlan)
    assert plan.step_policy == "coarse-only"
    assert plan.z0_ohm == pytest.approx(50.0)
    assert plan.port_order == ("P1", "P2", "P3", "P4")
    assert plan.frequencies_hz == pytest.approx((120.0e6, 125.0e6, 130.0e6))


def test_plan_frequency_sweep_coarse_plus_refined_merges_and_sorts_frequencies():
    plan = plan_frequency_sweep(
        start_hz=120.0e6,
        stop_hz=130.0e6,
        coarse_step_hz=5.0e6,
        refine_centers_hz=(127.0e6,),
        refine_half_span_hz=2.0e6,
        refined_step_hz=1.0e6,
        z0_ohm=75.0,
        port_order=("P1", "P2"),
    )

    assert plan.step_policy == "coarse+refined"
    assert plan.z0_ohm == pytest.approx(75.0)
    assert plan.port_order == ("P1", "P2")
    assert plan.frequencies_hz == pytest.approx(
        (
            120.0e6,
            125.0e6,
            126.0e6,
            127.0e6,
            128.0e6,
            129.0e6,
            130.0e6,
        )
    )


def test_plan_frequency_sweep_rejects_invalid_refined_config():
    with pytest.raises(ValueError, match="refined_step_hz must be positive"):
        plan_frequency_sweep(
            start_hz=120.0e6,
            stop_hz=130.0e6,
            coarse_step_hz=5.0e6,
            refine_centers_hz=(127.0e6,),
            refine_half_span_hz=1.0e6,
            refined_step_hz=0.0,
        )

    with pytest.raises(ValueError, match="refine_half_span_hz must be positive"):
        plan_frequency_sweep(
            start_hz=120.0e6,
            stop_hz=130.0e6,
            coarse_step_hz=5.0e6,
            refine_centers_hz=(127.0e6,),
            refine_half_span_hz=-1.0,
            refined_step_hz=1.0e6,
        )


# --- OPS-56: endpoint inclusion is decided relative to the step, not the carrier.
# Expected grids computed by hand: n = floor(span/step); the endpoint is appended
# iff span - n*step is not ~0 relative to the step. The pre-change module
# (default np.isclose, rtol 1e-5 of the carrier) returned 127.7000/127.7005 MHz,
# 127.7000/127.7010 MHz and a single 64 MHz point for the three cases below.
REPRODUCED_CASES = (
    ((127.7000e6, 127.7009e6, 500.0), (127.7000e6, 127.7005e6, 127.7009e6)),
    ((127.7000e6, 127.7012e6, 1.0e3), (127.7000e6, 127.7010e6, 127.7012e6)),
    ((64.0000e6, 64.0005e6, 1.0e3), (64.0000e6, 64.0005e6)),
)


@pytest.mark.parametrize("args, expected", REPRODUCED_CASES)
def test_uniform_grid_keeps_endpoint_at_khz_steps_on_mhz_carriers(args, expected):
    start, stop, step = args
    plan = plan_frequency_sweep(start_hz=start, stop_hz=stop, coarse_step_hz=step)
    grid = plan.frequencies_hz
    assert len(grid) == len(expected)
    # Absolute tolerance far below the step (1e-6 Hz), never carrier-relative.
    assert grid == pytest.approx(expected, rel=0.0, abs=1e-6)
    assert grid[-1] == stop  # endpoint present, bitwise


def test_uniform_grid_nondivisible_span_appends_stop():
    plan = plan_frequency_sweep(start_hz=0.0, stop_hz=10.0, coarse_step_hz=3.0)
    assert plan.frequencies_hz == (0.0, 3.0, 6.0, 9.0, 10.0)


def test_uniform_grid_float_noise_on_grid_endpoint_is_bitwise_stop():
    # 0.3 - 0.1 = 0.19999999999999998 and /0.1 = 1.9999999999999998: still 3 points.
    plan = plan_frequency_sweep(start_hz=0.1, stop_hz=0.3, coarse_step_hz=0.1)
    assert len(plan.frequencies_hz) == 3
    assert plan.frequencies_hz[-1] == 0.3  # bitwise, not approx
    assert plan.frequencies_hz == pytest.approx((0.1, 0.2, 0.3), rel=0.0, abs=1e-15)


def test_uniform_grid_on_grid_endpoint_bitwise_at_mhz():
    plan = plan_frequency_sweep(start_hz=63.0e6, stop_hz=65.0e6, coarse_step_hz=0.1e6)
    assert len(plan.frequencies_hz) == 21
    assert plan.frequencies_hz[-1] == 65.0e6


def test_coarse_refined_near_duplicate_merges_to_one_point():
    delta = 6.0e-8  # a few ulp at 125 MHz; << 1e-9 * refined step (1e-3 Hz)
    center = 127.0e6 + delta
    assert center - 2.0e6 != 125.0e6  # precondition: not bitwise equal
    plan = plan_frequency_sweep(
        start_hz=120.0e6,
        stop_hz=130.0e6,
        coarse_step_hz=5.0e6,
        refine_centers_hz=(center,),
        refine_half_span_hz=2.0e6,
        refined_step_hz=1.0e6,
    )
    grid = plan.frequencies_hz
    # np.unique would keep both 125e6 and 125e6+delta (8 points).
    assert len(grid) == 7
    assert grid[1] == 125.0e6  # first of the merged pair is kept
    assert grid == pytest.approx(
        (120.0e6, 125.0e6, 126.0e6, 127.0e6, 128.0e6, 129.0e6, 130.0e6), rel=0.0, abs=1e-6
    )


@pytest.mark.parametrize("field", ["start_hz", "stop_hz", "coarse_step_hz"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_uniform_grid_rejects_non_finite_inputs(field, bad):
    kwargs = {"start_hz": 120.0e6, "stop_hz": 130.0e6, "coarse_step_hz": 5.0e6}
    kwargs[field] = bad
    with pytest.raises(ValueError, match="must be finite"):
        plan_frequency_sweep(**kwargs)


# --- PORT-6: the planner equals the integer-lattice closed form, bitwise. ---
# On integer-Hz inputs every planner value start + k*step is an integer
# < 2**53, so float64 carries it exactly and the closed form (built in Python
# ints, never floats) must be matched bit for bit, not within a tolerance.

PORT6_SEED = 20260929
PORT6_N_CASES = 500
PORT6_COARSE_STEPS = (100, 500, 1_000, 10_000, 100_000, 1_000_000)
PORT6_REFINED_STEPS = (100, 250, 1_000)


def _int_lattice(start: int, stop: int, step: int) -> list[int]:
    """range(start, stop + 1, step) plus stop when stop is off the lattice."""
    points = list(range(start, stop + 1, step))
    if points[-1] != stop:
        points.append(stop)
    return points


def _int_reference(case: dict) -> list[int]:
    freqs = _int_lattice(case["start"], case["stop"], case["step"])
    for center in case["centers"]:
        lo = max(case["start"], center - case["half"])
        hi = min(case["stop"], center + case["half"])
        if hi < lo:
            continue
        freqs.extend(_int_lattice(lo, hi, case["rstep"]))
    return sorted(set(freqs))


def _port6_family() -> list[dict]:
    rng = np.random.default_rng(PORT6_SEED)
    cases = []
    for _ in range(PORT6_N_CASES):
        start = int(rng.integers(1_000_000, 300_000_001))
        step = int(rng.choice(PORT6_COARSE_STEPS))
        n_intervals = int(rng.integers(0, 61))
        remainder = 0 if rng.random() < 0.5 else int(rng.integers(1, step))
        stop = start + n_intervals * step + remainder
        n_centers = int(rng.integers(0, 4))
        centers: list[int] = []
        half = rstep = None
        if n_centers:
            allowed = [s for s in PORT6_REFINED_STEPS if s <= step]
            rstep = int(rng.choice(allowed))
            half = int(rng.integers(1, 10 * rstep + 1))
            # Centres may sit outside [start, stop]: exercises clipping and the
            # fully-outside (skipped) window.
            centers = [int(rng.integers(start - 2 * half, stop + 2 * half + 1)) for _ in range(n_centers)]
        z0 = float(rng.uniform(1.0, 200.0))
        port_order = tuple(f"P{i + 1}" for i in range(int(rng.integers(0, 5))))
        cases.append(
            dict(start=start, stop=stop, step=step, centers=centers, half=half, rstep=rstep,
                 z0=z0, port_order=port_order)
        )
    return cases


def _plan_kwargs(case: dict) -> dict:
    kwargs = dict(
        start_hz=float(case["start"]),
        stop_hz=float(case["stop"]),
        coarse_step_hz=float(case["step"]),
        z0_ohm=case["z0"],
        port_order=case["port_order"],
    )
    if case["centers"]:
        kwargs.update(
            refine_centers_hz=tuple(float(c) for c in case["centers"]),
            refine_half_span_hz=float(case["half"]),
            refined_step_hz=float(case["rstep"]),
        )
    return kwargs


def test_port6_planner_equals_integer_lattice_bitwise():
    cases = _port6_family()
    n_refined = sum(1 for c in cases if c["centers"])
    n_offgrid = sum(1 for c in cases if (c["stop"] - c["start"]) % c["step"])
    failures = []
    for idx, case in enumerate(cases):
        ref = tuple(map(float, _int_reference(case)))
        plan = plan_frequency_sweep(**_plan_kwargs(case))
        expected_policy = "coarse+refined" if case["centers"] else "coarse-only"
        ok = (
            plan.frequencies_hz == ref
            and plan.z0_ohm == case["z0"]
            and plan.port_order == case["port_order"]
            and plan.step_policy == expected_policy
        )
        if not ok:
            failures.append((idx, case, ref, plan))
    print(
        f"PORT6_FAMILY cases={len(cases)} refined={n_refined} offgrid_stop={n_offgrid} "
        f"failures={len(failures)}"
    )
    for idx, case, ref, plan in failures[:10]:
        print(f"PORT6_FAIL idx={idx} inputs={case} expected={ref} got={plan}")
    assert not failures


# Negative control: c1d167d's planner, transcribed verbatim from
# `git show c1d167d:src/fem_em_solver/ports/sweep.py` (only the two function
# names are prefixed). Prior record: 20260921T112323Z_OPS-56.log:34-37 (3/3 fail).
def _c1d167d_uniform_frequency_grid(start_hz: float, stop_hz: float, step_hz: float) -> np.ndarray:
    if step_hz <= 0.0:
        raise ValueError("step_hz must be positive")
    if stop_hz < start_hz:
        raise ValueError("stop_hz must be >= start_hz")

    span_hz = stop_hz - start_hz
    n_intervals = int(np.floor(span_hz / step_hz + 1e-12))

    frequencies = np.array([start_hz + idx * step_hz for idx in range(n_intervals + 1)], dtype=float)
    if not np.isclose(frequencies[-1], stop_hz):
        frequencies = np.append(frequencies, float(stop_hz))

    return frequencies


def _c1d167d_plan_frequency_sweep(
    *,
    start_hz: float,
    stop_hz: float,
    coarse_step_hz: float,
    refine_centers_hz: Sequence[float] | None = None,
    refine_half_span_hz: float | None = None,
    refined_step_hz: float | None = None,
    z0_ohm: float = 50.0,
    port_order: Sequence[str] = (),
) -> FrequencySweepPlan:
    """Plan a coarse-only or coarse+refined frequency sweep with deterministic grids."""
    if z0_ohm <= 0.0:
        raise ValueError("z0_ohm must be positive")

    coarse_grid = _c1d167d_uniform_frequency_grid(start_hz, stop_hz, coarse_step_hz)
    all_freqs: list[float] = list(coarse_grid)

    use_refined = bool(refine_centers_hz)
    if use_refined:
        if refine_half_span_hz is None or refine_half_span_hz <= 0.0:
            raise ValueError("refine_half_span_hz must be positive when refine_centers_hz are provided")
        if refined_step_hz is None or refined_step_hz <= 0.0:
            raise ValueError("refined_step_hz must be positive when refine_centers_hz are provided")

        for center_hz in refine_centers_hz or ():
            local_start = max(start_hz, float(center_hz) - refine_half_span_hz)
            local_stop = min(stop_hz, float(center_hz) + refine_half_span_hz)
            if local_stop < local_start:
                continue
            local_grid = _c1d167d_uniform_frequency_grid(local_start, local_stop, refined_step_hz)
            all_freqs.extend(local_grid.tolist())

    frequencies_hz = tuple(float(value) for value in np.unique(np.array(all_freqs, dtype=float)))
    step_policy = "coarse+refined" if use_refined else "coarse-only"

    return FrequencySweepPlan(
        frequencies_hz=frequencies_hz,
        step_policy=step_policy,
        z0_ohm=float(z0_ohm),
        port_order=tuple(port_order),
    )


C1D167D_LOGGED_WRONG_GRIDS = (
    ((127.7000e6, 127.7009e6, 500.0), (127.7000e6, 127.7005e6)),
    ((127.7000e6, 127.7012e6, 1.0e3), (127.7000e6, 127.7010e6)),
    ((64.0000e6, 64.0005e6, 1.0e3), (64.0000e6,)),
)


def test_port6_negative_control_c1d167d_reproduces_ops56_failures_and_misses_family():
    for (start, stop, step), wrong in C1D167D_LOGGED_WRONG_GRIDS:
        got = _c1d167d_plan_frequency_sweep(start_hz=start, stop_hz=stop, coarse_step_hz=step)
        print(f"PORT6_CONTROL case=({start!r}, {stop!r}, {step!r}) got={got.frequencies_hz}")
        assert got.frequencies_hz == wrong
    cases = _port6_family()
    n_fail = sum(
        1
        for case in cases
        if _c1d167d_plan_frequency_sweep(**_plan_kwargs(case)).frequencies_hz
        != tuple(map(float, _int_reference(case)))
    )
    print(f"PORT6_CONTROL_FAMILY c1d167d failures={n_fail}/{len(cases)} (predicted >= 1)")


# --- PORT-6 follow-up: the validity floor. ---------------------------------
# Below ~60 Hz steps the step-relative tolerance 1e-9*step falls under one ulp
# of a 100 MHz-class stop; the planner floors it at 4 ulp of the carrier, pins
# the merged grid's end to stop, and raises under 16 ulp.

PORT6_PROBE_CARRIERS = (63.8634e6 + 0.37, 127.7324e6 + 0.13, 297.2e6 + 0.71)


def _regime_probe(planner):
    """Non-integer carriers, 1 Hz / 0.1 Hz steps: one record per case.

    Duplicates are gaps < step/2; missing points are gaps > 1.5*step or a count
    short of n + 1.
    """
    records = []
    for carrier in PORT6_PROBE_CARRIERS:
        for step in (1.0, 0.1):
            for n, stop_form in itertools.product((7, 10, 37, 60), ("float-sum", "decimal")):
                start = carrier
                if stop_form == "float-sum":
                    stop = start + n * step
                else:  # stop as a user would type it: the exact decimal, rounded once
                    stop = float(Decimal(repr(start)) + n * Decimal(repr(step)))
                grid = planner(start_hz=start, stop_hz=stop, coarse_step_hz=step).frequencies_hz
                gaps = np.diff(np.asarray(grid))
                n_dup = int(np.sum(gaps < 0.5 * step))
                n_miss = int(np.sum(gaps > 1.5 * step)) + max(0, (n + 1) - (len(grid) - n_dup))
                records.append(
                    dict(carrier=carrier, step=step, n=n, stop_form=stop_form, start=start,
                         stop=stop, grid=grid, dup=n_dup, miss=n_miss)
                )
    return records


def _print_probe(label, records):
    for r in records:
        if r["dup"] or r["miss"]:
            gaps = np.diff(np.asarray(r["grid"]))
            print(
                f"{label} carrier={r['carrier']!r} step={r['step']} n={r['n']} stop={r['stop_form']} "
                f"len={len(r['grid'])} expected={r['n'] + 1} duplicates={r['dup']} missing={r['miss']} "
                f"min_gap={gaps.min()!r} last_two={r['grid'][-2:]}"
            )
    print(
        f"{label}_SUMMARY cases={len(records)} duplicates={sum(r['dup'] for r in records)} "
        f"missing={sum(r['miss'] for r in records)}"
    )


def test_port6_regime_probe_sub_60hz_steps():
    """(a) The 48-case probe, asserted: no duplicate, no missing point, exact ends."""
    records = _regime_probe(plan_frequency_sweep)
    _print_probe("PORT6_REGIME", records)
    assert len(records) == 48
    assert sum(r["dup"] for r in records) == 0
    assert sum(r["miss"] for r in records) == 0
    for r in records:
        assert len(r["grid"]) == r["n"] + 1, r
        assert r["grid"][0] == r["start"] and r["grid"][-1] == r["stop"], r


# (b) A sub-integer family against a Decimal reference: no float arithmetic is
# shared with the planner. random.Random, randrange / choice only, so the draw
# sequence is the same on every platform.
PORT6_FLOOR_SEED = 20260930
PORT6_FLOOR_N_CASES = 2000
PORT6_FLOOR_STEPS = tuple(
    Decimal(t) for t in ("0.01", "0.1", "0.25", "0.5", "1", "2.5", "10", "25", "50")
)
PORT6_FLOOR_ULP_BAND = 2


def _floor_family() -> list[dict]:
    rng = random.Random(PORT6_FLOOR_SEED)
    cases = []
    for _ in range(PORT6_FLOOR_N_CASES):
        carrier = rng.randrange(1_000_000, 500_000_000)
        cents = rng.randrange(100)
        start = Decimal(carrier) + Decimal(cents) / Decimal(100)
        step = rng.choice(PORT6_FLOOR_STEPS)
        n = rng.randrange(1, 201)
        kind = rng.choice(("float-sum", "decimal", "off-grid"))
        stop_dec = start + n * step
        if kind == "off-grid":
            stop_dec += rng.choice((Decimal("0.25"), Decimal("0.5"), Decimal("0.75"))) * step
        n_centers = rng.randrange(3)
        centers: list[Decimal] = []
        rstep = half = None
        if n_centers:
            rstep = rng.choice([t for t in PORT6_FLOOR_STEPS if t <= step])
            half = rstep * rng.randrange(1, 21)
            for _c in range(n_centers):
                j = rng.randrange(-2, n + 3)
                if rng.randrange(2):
                    offset = rstep * rng.randrange(0, 4)
                else:
                    offset = rng.choice((Decimal("0.003"), Decimal("0.017")))
                centers.append(start + j * step + offset)
        if kind == "float-sum":
            stop = float(start) + n * float(step)
        else:
            stop = float(stop_dec)
        cases.append(dict(start=start, step=step, n=n, kind=kind, stop_dec=stop_dec, stop=stop,
                          centers=centers, rstep=rstep, half=half))
    return cases


def _dec_lattice(lo: Decimal, hi: Decimal, step: Decimal) -> list[Decimal]:
    """lo + k*step <= hi, plus hi itself."""
    count = int((hi - lo) // step)
    return [lo + k * step for k in range(count + 1)] + [hi]


def _floor_reference(case: dict) -> list[Decimal]:
    points = _dec_lattice(case["start"], case["stop_dec"], case["step"])
    for center in case["centers"]:
        lo = max(case["start"], center - case["half"])
        hi = min(case["stop_dec"], center + case["half"])
        if hi < lo:
            continue
        points.extend(_dec_lattice(lo, hi, case["rstep"]))
    return sorted(set(points))


def _floor_kwargs(case: dict) -> dict:
    kwargs = dict(start_hz=float(case["start"]), stop_hz=case["stop"],
                  coarse_step_hz=float(case["step"]))
    if case["centers"]:
        kwargs.update(
            refine_centers_hz=tuple(float(c) for c in case["centers"]),
            refine_half_span_hz=float(case["half"]),
            refined_step_hz=float(case["rstep"]),
        )
    return kwargs


def _floor_check(planner, case: dict):
    """(ok, reason, worst deviation in ulp(stop)) of one case against the reference."""
    ref = _floor_reference(case)
    grid = planner(**_floor_kwargs(case)).frequencies_hz
    if len(grid) != len(ref):
        return False, f"len {len(grid)} != {len(ref)}", None
    if any(b <= a for a, b in zip(grid, grid[1:])):
        return False, "not strictly increasing", None
    if grid[0] != float(case["start"]) or grid[-1] != case["stop"]:
        return False, "ends not bitwise", None
    ulp = Fraction(math.ulp(case["stop"]))
    worst = max(abs(Fraction(g) - Fraction(r)) for g, r in zip(grid, ref)) / ulp
    if worst > PORT6_FLOOR_ULP_BAND:
        return False, f"deviation {float(worst):.3f} ulp", float(worst)
    return True, "", float(worst)


def test_port6_floor_family_equals_decimal_reference():
    """(b) 2 000 sub-integer cases: same length, increasing, exact ends, <= 2 ulp."""
    cases = _floor_family()
    failures = []
    worst = 0.0
    for idx, case in enumerate(cases):
        ok, reason, dev = _floor_check(plan_frequency_sweep, case)
        if dev is not None:
            worst = max(worst, dev)
        if not ok:
            failures.append((idx, reason, case))
    kinds = {k: sum(1 for c in cases if c["kind"] == k) for k in ("float-sum", "decimal", "off-grid")}
    print(
        f"PORT6_FLOOR_FAMILY cases={len(cases)} kinds={kinds} "
        f"refined={sum(1 for c in cases if c['centers'])} failures={len(failures)} "
        f"worst_deviation_ulp={worst:.3f} (predicted 1.456; band {PORT6_FLOOR_ULP_BAND})"
    )
    for idx, reason, case in failures[:10]:
        got = plan_frequency_sweep(**_floor_kwargs(case)).frequencies_hz
        print(f"PORT6_FLOOR_FAIL idx={idx} reason={reason} inputs={case} "
              f"expected={_floor_reference(case)} got={got}")
    assert not failures


def test_port6_floor_raises_under_16_ulp_and_is_exact_at_16():
    """(c) step 15 ulp raises; step 16 ulp gives 61 points with every gap 16 ulp."""
    start = 297200000.71
    u = float(np.spacing(start))
    with pytest.raises(ValueError, match="below the float64 resolution"):
        plan_frequency_sweep(start_hz=start, stop_hz=start + 60 * 15 * u, coarse_step_hz=15 * u)
    step = 16 * u
    grid = plan_frequency_sweep(
        start_hz=start, stop_hz=start + 60 * step, coarse_step_hz=step
    ).frequencies_hz
    gaps = np.diff(np.asarray(grid))
    print(f"PORT6_FLOOR_ULP u={u!r} step16={step!r} len={len(grid)} "
          f"gaps_all_16u={bool(np.all(gaps == step))}")
    assert len(grid) == 61
    assert np.all(gaps == step)


# Negative control: the pre-fix planner, transcribed verbatim from
# `git show 757e64ff:src/fem_em_solver/ports/sweep.py` (only the names are
# prefixed). Prior record: 20260930T020839Z_PORT-6-W1.log:69-71 (2 duplicates
# in 48: carrier 297200000.71, step 0.1, n = 7 and 37, decimal stop).
_757e64ff_STEP_RELATIVE_TOL = 1e-9


def _757e64ff_require_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return value


def _757e64ff_uniform_frequency_grid(start_hz: float, stop_hz: float, step_hz: float) -> np.ndarray:
    start_hz = _757e64ff_require_finite("start_hz", start_hz)
    stop_hz = _757e64ff_require_finite("stop_hz", stop_hz)
    step_hz = _757e64ff_require_finite("step_hz", step_hz)
    if step_hz <= 0.0:
        raise ValueError("step_hz must be positive")
    if stop_hz < start_hz:
        raise ValueError("stop_hz must be >= start_hz")

    span_hz = stop_hz - start_hz
    n_intervals = int(np.floor(span_hz / step_hz + _757e64ff_STEP_RELATIVE_TOL))

    frequencies = start_hz + step_hz * np.arange(n_intervals + 1, dtype=float)
    if abs(span_hz - n_intervals * step_hz) <= _757e64ff_STEP_RELATIVE_TOL * step_hz:
        # Endpoint is on the grid: pin it to the requested value exactly.
        frequencies[-1] = stop_hz
    else:
        frequencies = np.append(frequencies, stop_hz)

    return frequencies


def _757e64ff_merge_frequency_grids(values: Sequence[float], min_step_hz: float) -> tuple[float, ...]:
    """Sort and de-duplicate with a tolerance relative to the finest step in play."""
    tol_hz = _757e64ff_STEP_RELATIVE_TOL * min_step_hz
    merged: list[float] = []
    for value in np.sort(np.asarray(values, dtype=float)):
        if merged and value - merged[-1] <= tol_hz:
            continue
        merged.append(float(value))
    return tuple(merged)


def _757e64ff_plan_frequency_sweep(
    *,
    start_hz: float,
    stop_hz: float,
    coarse_step_hz: float,
    refine_centers_hz: Sequence[float] | None = None,
    refine_half_span_hz: float | None = None,
    refined_step_hz: float | None = None,
    z0_ohm: float = 50.0,
    port_order: Sequence[str] = (),
) -> FrequencySweepPlan:
    """Plan a coarse-only or coarse+refined frequency sweep with deterministic grids."""
    if z0_ohm <= 0.0:
        raise ValueError("z0_ohm must be positive")

    coarse_grid = _757e64ff_uniform_frequency_grid(start_hz, stop_hz, coarse_step_hz)
    all_freqs: list[float] = list(coarse_grid)
    min_step_hz = float(coarse_step_hz)

    use_refined = bool(refine_centers_hz)
    if use_refined:
        if refine_half_span_hz is None or refine_half_span_hz <= 0.0:
            raise ValueError("refine_half_span_hz must be positive when refine_centers_hz are provided")
        if refined_step_hz is None or refined_step_hz <= 0.0:
            raise ValueError("refined_step_hz must be positive when refine_centers_hz are provided")
        refine_half_span_hz = _757e64ff_require_finite("refine_half_span_hz", refine_half_span_hz)
        refined_step_hz = _757e64ff_require_finite("refined_step_hz", refined_step_hz)
        min_step_hz = min(min_step_hz, refined_step_hz)

        for center_hz in refine_centers_hz or ():
            center_hz = _757e64ff_require_finite("refine_centers_hz entry", center_hz)
            local_start = max(start_hz, center_hz - refine_half_span_hz)
            local_stop = min(stop_hz, center_hz + refine_half_span_hz)
            if local_stop < local_start:
                continue
            local_grid = _757e64ff_uniform_frequency_grid(local_start, local_stop, refined_step_hz)
            all_freqs.extend(local_grid.tolist())

    frequencies_hz = _757e64ff_merge_frequency_grids(all_freqs, min_step_hz)
    step_policy = "coarse+refined" if use_refined else "coarse-only"

    return FrequencySweepPlan(
        frequencies_hz=frequencies_hz,
        step_policy=step_policy,
        z0_ohm=float(z0_ohm),
        port_order=tuple(port_order),
    )


def test_port6_floor_negative_control_757e64ff_reproduces_logged_duplicates():
    records = _regime_probe(_757e64ff_plan_frequency_sweep)
    _print_probe("PORT6_FLOOR_CONTROL", records)
    bad = [(r["carrier"], r["step"], r["n"], r["stop_form"], r["dup"], r["miss"])
           for r in records if r["dup"] or r["miss"]]
    assert bad == [
        (297.2e6 + 0.71, 0.1, 7, "decimal", 1, 0),
        (297.2e6 + 0.71, 0.1, 37, "decimal", 1, 0),
    ]
    # Printed only (predicted 315/2000 wrong lengths; host arithmetic, not a record).
    cases = _floor_family()
    n_len = n_fail = 0
    for case in cases:
        ok, reason, _dev = _floor_check(_757e64ff_plan_frequency_sweep, case)
        n_fail += not ok
        n_len += reason.startswith("len")
    print(f"PORT6_FLOOR_CONTROL_FAMILY 757e64ff wrong_lengths={n_len}/{len(cases)} "
          f"(predicted 315), any_failure={n_fail}")
