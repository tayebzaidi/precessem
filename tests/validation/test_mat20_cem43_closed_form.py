"""`MAT-20` step 1 — CEM43 thermal dose against its closed form (§9 item 318).

The operator ``thermal.dose.cem43`` integrates R^(43 − T) exactly on the
piecewise-linear interpolant of a sampled history (ruling 318-1,
docs/planning/chunks/MAT-20.md).  Anchors, all at the §7 row's 1e-10 relative:

(A1) three linear ramps straddling 43 °C, each at 2 / 7 / 1000 samples, against
     ``ramp_dose`` below — an independent two-piece closed form that does not
     call the operator (raw difference quotient, no expm1);
(A2) `MAT-9`'s closed-form step response at the ball centre, re-based to
     T∞ = T_init = 37 °C so that it crosses 43 °C, against ``scipy.integrate.quad``
     split at the ``brentq`` crossing; N, 2N, 4N … printed, asserted at the
     first N whose predicted (Richardson, rate 2) error is ≤ 1e-11.

Controls: (Rc) the planted R swap misses (A1) by > 1e-10 on every ramp
(asserted); (Sc) no split at 43 °C — printed only.  The plain rectangle sum
Σ R^(43−T_i)Δt_i is printed beside the operator, never asserted.

Scope: a post-processor against closed forms; no limit, threshold, tissue or
compliance claim, and no FEM field.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from fem_em_solver.thermal.dose import cem43
from fem_em_solver.utils.analytical import robin_sphere_transient_temperature

BAR = 1e-10
RAMPS = {  # name: (T_start, T_end, duration [s])
    "rise 41->45 / 600 s": (41.0, 45.0, 600.0),
    "shallow 42.5->43.5 / 60 s": (42.5, 43.5, 60.0),
    "fall 45->40 / 900 s": (45.0, 40.0, 900.0),
}
SAMPLINGS = (2, 7, 1000)

# MAT-9 constants (test_mat9_transient_conduction_sphere.py:41-45), re-based.
R, K, RHO_C, Q, H = 0.1, 0.6, 4.18e6, 2000.0, 10.0
T_BASE = 37.0
T_END = 36000.0
PHYS = dict(radius=R, k=K, rho_c=RHO_C, q=Q, h=H, t_inf=T_BASE, t_init=T_BASE)


def _side(ta, tb, tau, r):
    """∫₀^τ r^(43 − ta − b t) dt, b = (tb − ta)/τ, on one side of 43."""
    return (r ** (43.0 - ta) - r ** (43.0 - tb)) / ((tb - ta) / tau * math.log(r))


def ramp_dose(a, b_end, tau, r_above=0.5, r_below=0.25):
    """Analytic dose of the ramp a → b_end over τ, split at 43 °C."""
    r_of = lambda T: r_above if T >= 43.0 else r_below  # noqa: E731
    if (a - 43.0) * (b_end - 43.0) >= 0.0:
        return _side(a, b_end, tau, r_of(0.5 * (a + b_end)))
    t43 = (43.0 - a) / ((b_end - a) / tau)
    return _side(a, 43.0, t43, r_of(a)) + _side(43.0, b_end, tau - t43, r_of(b_end))


def _ramp(a, b_end, tau, n):
    t = np.linspace(0.0, tau, n)
    return t, a + (b_end - a) * t / tau


def _rect_sum(t, T, r_above=0.5, r_below=0.25):
    T0 = T[:-1]
    return float(np.sum(np.where(T0 >= 43.0, r_above, r_below) ** (43.0 - T0) * np.diff(t)))


def _rel(x, ref):
    return abs(x - ref) / abs(ref)


@pytest.mark.parametrize("name", list(RAMPS))
def test_a1_ramps_exact(name):
    a, b_end, tau = RAMPS[name]
    ref = ramp_dose(a, b_end, tau)
    for n in SAMPLINGS:
        t, T = _ramp(a, b_end, tau, n)
        got = cem43(t, T)
        err = _rel(got, ref)
        print(
            f"MAT20s1 A1 {name:28s} n={n:5d} ref={ref:.15e} op={got:.15e} "
            f"rel={err:.3e} rect_sum={_rect_sum(t, T):.10e}"
        )
        assert err <= BAR, (name, n, err)


def test_a1_two_d_columns():
    t = np.linspace(0.0, 600.0, 1000)
    ends = [(41.0, 45.0), (42.5, 43.5), (45.0, 40.0)]
    T = np.stack([a + (b - a) * t / 600.0 for a, b in ends], axis=1)
    got = cem43(t, T)
    assert got.shape == (3,)
    for j, (a, b) in enumerate(ends):
        err = _rel(got[j], ramp_dose(a, b, 600.0))
        print(f"MAT20s1 A1-2D col={j} rel={err:.3e}")
        assert err <= BAR


def test_rc_r_swap_misses():
    for name, (a, b_end, tau) in RAMPS.items():
        ref = ramp_dose(a, b_end, tau)
        for n in SAMPLINGS:
            t, T = _ramp(a, b_end, tau, n)
            miss = _rel(cem43(t, T, r_above=0.25, r_below=0.5), ref)
            print(f"MAT20s1 Rc {name:28s} n={n:5d} swap_miss={miss:.4f}")
            assert miss > BAR, (name, n, miss)


def test_sc_no_split_printed():
    for name, (a, b_end, tau) in RAMPS.items():
        ref = ramp_dose(a, b_end, tau)
        r = 0.5 if a >= 43.0 else 0.25  # one R (the start side) over the whole crossing segment
        nosplit = _side(a, b_end, tau, r)
        print(
            f"MAT20s1 Sc {name:28s} no_split_miss={_rel(nosplit, ref):.4f} (printed, not asserted)"
        )


def test_complex_nonzero_imag_raises():
    t = np.array([0.0, 1.0])
    with pytest.raises(ValueError):
        cem43(t, np.array([42.0 + 1e-3j, 44.0]))
    assert _rel(cem43(t, np.array([42.0 + 0j, 44.0 + 0j])), ramp_dose(42.0, 44.0, 1.0)) <= BAR


def _T_centre(t):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.empty_like(t)
    for s in range(0, t.size, 8192):
        out[s : s + 8192] = robin_sphere_transient_temperature(0.0, t[s : s + 8192], **PHYS)
    return out


def test_a2_mat9_step_response():
    T0, T1 = _T_centre([0.0, T_END])
    print(f"MAT20s1 A2 bracket T(0)={T0:.12f} T(T_END={T_END:g})={T1:.12f}")
    assert T0 < 43.0 < T1, "brentq bracket"
    t43 = brentq(
        lambda s: _T_centre(s)[0] - 43.0, 0.0, T_END, xtol=1e-12, rtol=4 * np.finfo(float).eps
    )

    def f(s):
        T = _T_centre(s)[0]
        return (0.5 if T >= 43.0 else 0.25) ** (43.0 - T)

    i1, e1 = quad(f, 0.0, t43, epsabs=1e-13, epsrel=1e-13, limit=2000)
    i2, e2 = quad(f, t43, T_END, epsabs=1e-13, epsrel=1e-13, limit=2000)
    ref = i1 + i2
    print(
        f"MAT20s1 A2 t43={t43:.10f} quad={ref:.15e} quad_err_est={e1 + e2:.3e} "
        f"(rel {(e1 + e2) / ref:.3e})"
    )

    k_max = 19
    t_fine = np.linspace(0.0, T_END, 2**k_max + 1)
    T_fine = _T_centre(t_fine)
    prev = prev2 = None
    asserted = None
    for k in range(10, k_max + 1):
        stride = 2 ** (k_max - k)
        t, T = t_fine[::stride], T_fine[::stride]
        d = cem43(t, T)
        err = _rel(d, ref)
        pred = abs(d - prev) / 3.0 / abs(d) if prev is not None else float("nan")
        rate = (
            math.log2(abs(prev - prev2) / abs(d - prev))
            if k >= 12 and d != prev and prev != prev2
            else float("nan")
        )
        print(
            f"MAT20s1 A2 N={2**k:7d} op={d:.15e} rel_vs_quad={err:.3e} "
            f"pred={pred:.3e} rate={rate:.3f} rect_sum={_rect_sum(t, T):.10e}"
        )
        if asserted is None and pred <= 1e-11:
            asserted = (2**k, err)
        prev2, prev = prev, d
        if asserted is not None:
            break
    assert asserted is not None, "no N up to 2^19 reached predicted error <= 1e-11"
    print(f"MAT20s1 A2 asserted at N={asserted[0]} rel={asserted[1]:.3e}")
    assert asserted[1] <= BAR, asserted
