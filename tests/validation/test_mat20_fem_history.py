"""`MAT-20` step 2, FEM clause — the six-minute operator on `MAT-17` step 3's history.

§9 item 370 (chain P19; ruling 319-3 option (a), ratified by the 2026-10-10
weekly).  `MAT-17` step 3's fixture (`MAT-15` step 1's `_run_heated`, F-small
coil-driven gel phantom, 10 MHz, k = 0, G = 1 W / P₁) is re-run at `-n 2` with
every step's volume-mean ⟨θ⟩ₙ kept, for the duty schedule and the ten
full-from-zero Heaviside runs H_m.

Asserted (bands imported): (H0) the per-step series obeys Duhamel,
max |⟨θ_duty⟩ₙ − Σ(−1)^m ⟨θ_Hm⟩ₙ| / max|⟨θ_duty⟩| ≤ ``L_BAND``; (FEM) the
operator's 360 s window mean of the duty series equals Σ(−1)^m of it on the
H_m series ≤ ``D_BAND``; (FEM-2d) the same through the 2-D call.
Negative control (Wc-FEM): a shifted by P/2 (``STEPS_PER_HALF`` samples)
against b misses by > ``D_BAND``.  Printed only: (Dc) m = 9 dropped, the
step-3 log's end-time records, t_max and the maximum six-minute means.

Scope: an operator identity on one FEM history.  Not: a limit value,
compliance, tissue, perfusion, an absolute temperature.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.thermal.dose import running_average
from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    N_MAIN,
    T_INF,
    T_STAR,
    _run_heated,
)
from tests.validation.test_mat17_duty_cycle_coil_driven import (
    D_BAND,
    DT,
    HALF,
    L_BAND,
    N_SWITCH,
    STEPS_PER_HALF,
    _chi_duty,
    _heaviside,
    _transient_s,
)

W = 360.0
OUT_JSON = "logs/mat20/mat17s3-series-n2.json"
STEPS = tuple(range(N_MAIN + 1))


def _series(res, work, f_mean, comm, v_sub):
    out = []
    for n in STEPS:
        work.x.array[:] = res.snapshots[n].x.array - T_INF
        work.x.scatter_forward()
        out.append(float(np.real(comm.allreduce(fem.assemble_scalar(f_mean), op=MPI.SUM))) / v_sub)
    res.snapshots.clear()
    return np.asarray(out)


def _fem_history(*, keep_duty_fields=False):
    """`hist()`'s body (rule (a) lift, EX-110): the fixture run, no JSON dump.

    ``keep_duty_fields`` additionally returns the fixture objects and the duty
    run's per-step P1 fields (theta = T - T_INF, copied before ``_series``
    clears the snapshots); the default return is unchanged.
    """
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    base = _run_heated(return_fields=True)
    t_base = time.perf_counter() - t0
    sub, ft, v_sub = base["sub"], base["ft"], base["v_sub"]
    q0 = base["q_h"][0]
    V1 = fem.functionspace(sub, ("Lagrange", 1))
    work = fem.Function(V1)
    f_mean = fem.form(work * ufl.dx(domain=sub))

    t1 = time.perf_counter()
    r_duty = _transient_s(sub, ft, q0, _chi_duty, snapshot_steps=STEPS)
    duty_fields = None
    if keep_duty_fields:
        duty_fields = [np.real(r_duty.snapshots[n].x.array).copy() - T_INF for n in STEPS]
    duty = _series(r_duty, work, f_mean, comm, v_sub)
    hs = [
        _series(
            _transient_s(sub, ft, q0, _heaviside(m), snapshot_steps=STEPS),
            work,
            f_mean,
            comm,
            v_sub,
        )
        for m in range(N_SWITCH)
    ]
    t_thermal = time.perf_counter() - t1
    times = np.arange(N_MAIN + 1) * DT
    _say(
        f"\nMAT20fem fixture: _run_heated {t_base:.1f} s; 11 transients + "
        f"{11 * len(STEPS)} reductions {t_thermal:.1f} s; DT = {DT:.6f} s, "
        f"W/DT = {W / DT:.4f} samples, ranks = {comm.size}"
    )
    out = {"t": times, "duty": duty, "H": hs}
    if keep_duty_fields:
        out.update(
            {
                "base": base,
                "sub": sub,
                "ft": ft,
                "V1": V1,
                "work": work,
                "f_mean": f_mean,
                "v_sub": v_sub,
                "q0": q0,
                "duty_fields": duty_fields,
            }
        )
    return out


@pytest.fixture(scope="module")
def hist():
    comm = MPI.COMM_WORLD
    out = _fem_history()
    times, duty, hs = out["t"], out["duty"], out["H"]
    if comm.rank == 0:
        os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
        with open(OUT_JSON, "w") as fh:
            json.dump(
                {
                    "times": times.tolist(),
                    "duty": duty.tolist(),
                    "H": [h.tolist() for h in hs],
                    "DT": DT,
                    "HALF": HALF,
                    "N_SWITCH": N_SWITCH,
                    "T_INF": T_INF,
                    "T_STAR": T_STAR,
                    "ranks": comm.size,
                },
                fh,
                indent=1,
            )
    comm.Barrier()
    return out


def _identities(t, duty, hs):
    """The test body's computations (rule (a) lift, EX-110); prints stay in the test."""
    signs = [(-1.0) ** m for m in range(N_SWITCH)]
    sup = sum(s * h for s, h in zip(signs, hs))
    h0 = float(np.max(np.abs(duty - sup)) / np.max(np.abs(duty)))
    a, a_mx, a_tmx = running_average(t, duty, window=W)
    ra = [running_average(t, h, window=W) for h in hs]
    b = sum(s * r[0] for s, r in zip(signs, ra))
    ok = ~np.isnan(a)
    scale = np.max(np.abs(a[ok]))
    res = np.abs(a[ok] - b[ok]) / scale
    fem_v = float(np.max(res))
    arr = np.column_stack([duty] + hs)
    m2, _, _ = running_average(t, arr, window=W)
    b2 = m2[:, 1:] @ np.asarray(signs)
    fem2d = float(np.max(np.abs(m2[ok, 0] - b2[ok])) / np.max(np.abs(m2[ok, 0])))
    fem2d_eq = bool(np.array_equal(m2[ok, 0], a[ok]))
    idx = np.nonzero(ok)[0]
    idx = idx[idx - STEPS_PER_HALF >= idx[0]]
    wc = float(np.max(np.abs(a[idx] - b[idx - STEPS_PER_HALF])) / scale)
    b9 = b - signs[9] * ra[9][0]
    dc = float(np.max(np.abs(a[ok] - b9[ok])) / scale)
    _, h0_mx, h0_tmx = ra[0]
    return {
        "signs": signs,
        "h0": h0,
        "a": a,
        "a_mx": a_mx,
        "a_tmx": a_tmx,
        "ra": ra,
        "b": b,
        "ok": ok,
        "scale": scale,
        "res": res,
        "fem": fem_v,
        "fem2d": fem2d,
        "fem2d_eq": fem2d_eq,
        "idx": idx,
        "wc": wc,
        "dc": dc,
        "h0_mx": h0_mx,
        "h0_tmx": h0_tmx,
    }


@complex_only
def test_mat20_fem_history(hist):
    t, duty, hs = hist["t"], hist["duty"], hist["H"]
    ids = _identities(t, duty, hs)

    # Records against the step-3 W1 log (:1960–1969, :1974).
    for m, h in enumerate(hs):
        _say(f"MAT20fem record H_{m} (t_m = {m * HALF:.4f} s): <theta>(t*) = {h[-1]:.12e} K")
    _say(
        f"MAT20fem record duty <theta> n=32 {duty[32]:.9f} K  n=36 {duty[36]:.9f} K  "
        f"n=40 {duty[40]:.9f} K  (step-3 log :1974: 0.466250330 / 0.584438780 / 0.568835514)"
    )

    # (H0) Duhamel on the per-step series.
    h0 = ids["h0"]
    _say(
        f"MAT20fem (H0) max_n |duty - sum (-1)^m H_m| / max|duty| = {h0:.3e} "
        f"(bound L_BAND {L_BAND:.0e}, predicted <= 5e-14)"
    )

    # (FEM) the operator column-wise.
    fem_v = ids["fem"]
    _say(
        f"MAT20fem (FEM) W = {W} s: max |a - b| / max|a| over {int(ids['ok'].sum())} windows = "
        f"{fem_v:.3e} (bound D_BAND {D_BAND:.0e}, predicted <= 5e-14)"
    )
    _say("MAT20fem (FEM) per-window residuals: " + " ".join(f"{x:.2e}" for x in ids["res"]))

    # (FEM-2d) through the 2-D call.
    fem2d = ids["fem2d"]
    _say(
        f"MAT20fem (FEM-2d) 2-D call: {fem2d:.3e} (bound D_BAND {D_BAND:.0e}); "
        f"2-D duty column == 1-D: {ids['fem2d_eq']}"
    )

    # (Wc-FEM) a at t_i against b at t_{i-4}.
    wc = ids["wc"]
    _say(
        f"MAT20fem control (Wc-FEM) shift P/2 = {STEPS_PER_HALF} samples: miss {wc:.3e} "
        f"over {ids['idx'].size} windows (asserted > {D_BAND:.0e}; predicted ~1e-2-1e-1)"
    )

    # (Dc) printed only: m = 9 dropped.
    _say(f"MAT20fem control (Dc) m = 9 dropped, printed: {ids['dc']:.6f}")

    # Records: six-minute maxima.
    _say(
        f"MAT20fem record duty max six-minute mean {ids['a_mx']:.9e} K at t_max = "
        f"{ids['a_tmx']:.4f} s; H_0 max {ids['h0_mx']:.9e} K at t_max = "
        f"{ids['h0_tmx']:.4f} s (t* = {T_STAR:.4f} s)"
    )

    assert h0 <= L_BAND
    assert fem_v <= D_BAND
    assert fem2d <= D_BAND
    assert wc > D_BAND
