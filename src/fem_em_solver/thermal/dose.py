"""Thermal-dose post-processors (Phase 8, `MAT-20`).

Step 1: :func:`cem43`, the cumulative equivalent minutes at 43 °C of a
temperature history (Sapareto & Dewey's CEM43 with R = 0.5 above and
R = 0.25 below the reference), evaluated **exactly on the piecewise-linear
interpolant** of the samples (ruling 318-1, `docs/planning/chunks/MAT-20.md`).

Between samples ``T`` is linear in ``t``; on such a segment

    ∫ R^(t_ref − T) dt = Δt · R^(t_ref − T₀) · expm1(x) / x,   x = −ΔT ln R,

which is the closed form ``(R^(t_ref−T₀) − R^(t_ref−T₁)) / (b ln R)``
(``b = ΔT/Δt``) written so that it does not cancel as ``ΔT → 0``.  A segment
that crosses ``t_ref`` is split at ``t₀ + (t_ref − T₀)/b`` and each part uses
the R of its own side; the integrand equals 1 on both sides at ``t_ref``.

The result is in the units of ``times`` (seconds in → equivalent seconds at
43 °C out; divide by 60 for CEM43 minutes).  No dose limit, damage threshold
or tissue model is implied.

Step 2: :func:`running_average`, the trailing window mean (default 360 s)
of a sampled series, exact on the same piecewise-linear interpolant.
"""

from __future__ import annotations

import numpy as np

__all__ = ["cem43", "running_average"]


def _phi(x):
    """``expm1(x)/x`` with the removable singularity at 0 filled by 1."""
    x = np.asarray(x, dtype=float)
    small = x == 0.0
    safe = np.where(small, 1.0, x)
    return np.where(small, 1.0, np.expm1(safe) / safe)


def _piece(t_a, t_b, dur, ln_r, t_ref):
    """Exact ∫ R^(t_ref − T) over a linear piece from ``t_a`` to ``t_b``."""
    return dur * np.exp(ln_r * (t_ref - t_a)) * _phi(-ln_r * (t_b - t_a))


def cem43(times, temps, *, r_above=0.5, r_below=0.25, t_ref=43.0):
    """Thermal dose ∫ R^(t_ref − T(t)) dt of a sampled history.

    Parameters
    ----------
    times : (n_t,) array, strictly increasing.
    temps : (n_t,) or (n_t, n_cells) array of temperatures [°C] — e.g. one
        DG0 value per cell per time step.  A complex array is accepted only
        if every imaginary part is exactly zero; otherwise ``ValueError``
        (a complex-mode history is never silently truncated to its real part).
    r_above, r_below : R on ``T ≥ t_ref`` and ``T < t_ref``.
    t_ref : reference temperature.

    Returns
    -------
    float for 1-D ``temps``; ``(n_cells,)`` array for 2-D.
    """
    t = np.asarray(times)
    T = np.asarray(temps)
    if np.iscomplexobj(t):
        raise ValueError("cem43: times must be real")
    if np.iscomplexobj(T):
        if np.any(T.imag != 0.0):
            raise ValueError("cem43: temperatures have a non-zero imaginary part")
        T = T.real  # exactly zero imaginary part: nothing is discarded
    t = t.astype(float)
    T = T.astype(float)
    if t.ndim != 1 or t.size < 2:
        raise ValueError("cem43: times must be 1-D with at least two samples")
    if T.ndim not in (1, 2) or T.shape[0] != t.size:
        raise ValueError("cem43: temps must be (n_t,) or (n_t, n_cells)")
    dt = np.diff(t)
    if np.any(dt <= 0.0):
        raise ValueError("cem43: times must be strictly increasing")
    if r_above <= 0.0 or r_below <= 0.0:
        raise ValueError("cem43: R must be positive")

    one_d = T.ndim == 1
    if one_d:
        T = T[:, None]
    dt = dt[:, None]
    T0, T1 = T[:-1], T[1:]
    ln_a, ln_b = np.log(r_above), np.log(r_below)

    cross = (T0 - t_ref) * (T1 - t_ref) < 0.0
    dT = T1 - T0
    frac = np.where(cross, (t_ref - T0) / np.where(cross, dT, 1.0), 1.0)
    # Piece 1: T0 → (t_ref if crossing else T1).  Its side is T0's when the
    # segment crosses, else the side of the segment mean (a segment touching
    # t_ref at one end lies wholly on the other end's side).
    end1 = np.where(cross, t_ref, T1)
    above1 = np.where(cross, T0 > t_ref, 0.5 * (T0 + T1) >= t_ref)
    ln1 = np.where(above1, ln_a, ln_b)
    dose = _piece(T0, end1, frac * dt, ln1, t_ref)
    # Piece 2 (crossing segments only): t_ref → T1, on T1's side.
    ln2 = np.where(T1 > t_ref, ln_a, ln_b)
    dose2 = _piece(np.full_like(T1, t_ref), T1, (1.0 - frac) * dt, ln2, t_ref)
    dose = dose + np.where(cross, dose2, 0.0)

    out = dose.sum(axis=0)
    return float(out[0]) if one_d else out


def running_average(times, values, *, window=360.0):
    """Trailing window mean ``(1/W) ∫_{t−W}^{t} v(s) ds`` of a sampled series.

    ``v`` is the **piecewise-linear interpolant** of the samples and the
    integral is exact on it, including the fractional segment at the window's
    start: ``W`` need not be a whole number of samples (`MAT-20` step 2).
    The mean is formed at every sample time ``tᵢ`` with ``tᵢ − W ≥ t₀``;
    earlier entries are NaN.

    Parameters
    ----------
    times : (n_t,) strictly increasing.
    values : (n_t,) or (n_t, n_cells) real array (same complex rule as
        :func:`cem43`).
    window : W > 0, in the units of ``times``; must not exceed the span.

    Returns
    -------
    (means, max, t_max): ``means`` has the shape of ``values``; ``max`` and
    ``t_max`` are floats for 1-D input and ``(n_cells,)`` arrays for 2-D.
    No exposure limit or compliance statement is implied.
    """
    t = np.asarray(times)
    v = np.asarray(values)
    if np.iscomplexobj(t):
        raise ValueError("running_average: times must be real")
    if np.iscomplexobj(v):
        if np.any(v.imag != 0.0):
            raise ValueError("running_average: values have a non-zero imaginary part")
        v = v.real
    t = t.astype(float)
    v = v.astype(float)
    if t.ndim != 1 or t.size < 2:
        raise ValueError("running_average: times must be 1-D with at least two samples")
    if v.ndim not in (1, 2) or v.shape[0] != t.size:
        raise ValueError("running_average: values must be (n_t,) or (n_t, n_cells)")
    dt = np.diff(t)
    if np.any(dt <= 0.0):
        raise ValueError("running_average: times must be strictly increasing")
    W = float(window)
    if not W > 0.0 or W > t[-1] - t[0]:
        raise ValueError("running_average: need 0 < window <= times[-1] - times[0]")

    one_d = v.ndim == 1
    if one_d:
        v = v[:, None]
    # Exact cumulative integral of the interpolant at the samples.
    F = np.zeros_like(v)
    F[1:] = np.cumsum(0.5 * (v[1:] + v[:-1]) * dt[:, None], axis=0)

    valid = t - W >= t[0]
    idx = np.nonzero(valid)[0]
    s = t[idx] - W  # window starts
    # Segment j with t_j <= s < t_{j+1} (s == t_last cannot occur: W > 0).
    j = np.clip(np.searchsorted(t, s, side="right") - 1, 0, t.size - 2)
    a = (s - t[j]) / dt[j]  # fraction into segment j
    v_s = v[j] + a[:, None] * (v[j + 1] - v[j])
    F_s = F[j] + (s - t[j])[:, None] * 0.5 * (v[j] + v_s)

    means = np.full_like(v, np.nan)
    means[idx] = (F[idx] - F_s) / W
    k = np.argmax(means[idx], axis=0)
    mx = means[idx][k, np.arange(v.shape[1])]
    t_mx = t[idx][k]
    if one_d:
        return means[:, 0], float(mx[0]), float(t_mx[0])
    return means, mx, t_mx
