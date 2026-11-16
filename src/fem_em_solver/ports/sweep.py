"""Frequency sweep planning helpers (chunk D4)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class FrequencySweepPlan:
    """Deterministic frequency sweep plan + metadata for traceability."""

    frequencies_hz: tuple[float, ...]
    step_policy: str
    z0_ohm: float
    port_order: tuple[str, ...]


# Grid tolerances are relative to the STEP, never to the carrier (OPS-56):
# a default ``np.isclose`` (rtol 1e-5) at 128 MHz treated anything within
# ~1.3 kHz of ``stop_hz`` as already present, so kHz refinement steps lost
# their endpoint (``(64.0e6, 64.0005e6, 1e3)`` returned a single point).
_STEP_RELATIVE_TOL = 1e-9
# ...floored at a few ulp of the carrier (PORT-6 follow-up): under ~60 Hz
# steps ``1e-9 * step`` is below one ulp of a 100 MHz-class ``stop``, so a
# decimal-typed on-grid ``stop`` one ulp from ``start + n * step`` was kept as
# a near-duplicate (2 of 48 probe cases, 20260930T020839Z_PORT-6-W1.log:69-71).
_CARRIER_ULP_TOL = 4.0
# A step under this many ulp of the carrier is not a float64 grid at all.
_MIN_STEP_ULPS = 16.0


def _grid_tol_hz(step_hz: float, at_hz: float) -> float:
    return max(_STEP_RELATIVE_TOL * step_hz, _CARRIER_ULP_TOL * float(np.spacing(abs(at_hz))))


def _require_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return value


def _uniform_frequency_grid(start_hz: float, stop_hz: float, step_hz: float) -> np.ndarray:
    start_hz = _require_finite("start_hz", start_hz)
    stop_hz = _require_finite("stop_hz", stop_hz)
    step_hz = _require_finite("step_hz", step_hz)
    if step_hz <= 0.0:
        raise ValueError("step_hz must be positive")
    if stop_hz < start_hz:
        raise ValueError("stop_hz must be >= start_hz")
    carrier_hz = max(abs(start_hz), abs(stop_hz))
    min_step_hz = _MIN_STEP_ULPS * float(np.spacing(carrier_hz))
    if step_hz < min_step_hz:
        raise ValueError(
            f"step_hz={step_hz!r} is below the float64 resolution of the grid: "
            f"{_MIN_STEP_ULPS:g} ulp at {carrier_hz!r} Hz is {min_step_hz!r} Hz"
        )

    span_hz = stop_hz - start_hz
    n_intervals = int(np.floor(span_hz / step_hz + _STEP_RELATIVE_TOL))

    frequencies = start_hz + step_hz * np.arange(n_intervals + 1, dtype=float)
    if abs(span_hz - n_intervals * step_hz) <= _grid_tol_hz(step_hz, carrier_hz):
        # Endpoint is on the grid: pin it to the requested value exactly.
        frequencies[-1] = stop_hz
    else:
        frequencies = np.append(frequencies, stop_hz)

    return frequencies


def _merge_frequency_grids(values: Sequence[float], min_step_hz: float) -> tuple[float, ...]:
    """Sort and de-duplicate with a tolerance relative to the finest step in play.

    A cluster keeps its first member; the last cluster keeps the largest value,
    so the grid still ends bitwise on ``stop_hz``.
    """
    ordered = np.sort(np.asarray(values, dtype=float))
    merged: list[float] = []
    for value in ordered:
        if merged and value - merged[-1] <= _grid_tol_hz(min_step_hz, value):
            continue
        merged.append(float(value))
    if merged:
        merged[-1] = float(ordered[-1])
    return tuple(merged)


def plan_frequency_sweep(
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

    coarse_grid = _uniform_frequency_grid(start_hz, stop_hz, coarse_step_hz)
    all_freqs: list[float] = list(coarse_grid)
    min_step_hz = float(coarse_step_hz)

    use_refined = bool(refine_centers_hz)
    if use_refined:
        if refine_half_span_hz is None or refine_half_span_hz <= 0.0:
            raise ValueError("refine_half_span_hz must be positive when refine_centers_hz are provided")
        if refined_step_hz is None or refined_step_hz <= 0.0:
            raise ValueError("refined_step_hz must be positive when refine_centers_hz are provided")
        refine_half_span_hz = _require_finite("refine_half_span_hz", refine_half_span_hz)
        refined_step_hz = _require_finite("refined_step_hz", refined_step_hz)
        min_step_hz = min(min_step_hz, refined_step_hz)

        for center_hz in refine_centers_hz or ():
            center_hz = _require_finite("refine_centers_hz entry", center_hz)
            local_start = max(start_hz, center_hz - refine_half_span_hz)
            local_stop = min(stop_hz, center_hz + refine_half_span_hz)
            if local_stop < local_start:
                continue
            local_grid = _uniform_frequency_grid(local_start, local_stop, refined_step_hz)
            all_freqs.extend(local_grid.tolist())

    frequencies_hz = _merge_frequency_grids(all_freqs, min_step_hz)
    step_policy = "coarse+refined" if use_refined else "coarse-only"

    return FrequencySweepPlan(
        frequencies_hz=frequencies_hz,
        step_policy=step_policy,
        z0_ohm=float(z0_ohm),
        port_order=tuple(port_order),
    )
